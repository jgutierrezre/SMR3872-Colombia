# PyTorch Lightning

[PyTorch Lightning](https://www.pytorchlightning.ai) wraps PyTorch to provide easy, distributed training done in your choice of numerical precision for the matrix multiples and other operations. To convert from PyTorch to PyTorch Lightning one simply needs to:

+ restructure the code by moving the network definition, optimizer and other code to a subclass of `pl.LightningModule`  
+ remove `.cuda()` and `.to()` calls since Lightning code is hardware agnostic  

Once these changes have been made one can simply choose how many nodes or GPUs to use and Lightning will take care of the rest. One can also use different numerical precisions (fp16, bf16). There is tensorboard support.

## Requirements

To install  lightning we can use pip:

$ pip install pytorch-lightning


See the [Trainer API](https://pytorch-lightning.readthedocs.io/en/latest/common/trainer.html#trainer-class-api).

## Single-GPU Example

Below is an example PL script:

```python
import torch
from torch import nn
from torch.nn import functional as F
from torch.utils.data import DataLoader
from torch.utils.data import random_split
from torchvision.datasets import MNIST
from torchvision import transforms
import pytorch_lightning as pl
import os

class LitAutoEncoder(pl.LightningModule):
    def __init__(self):
        super().__init__()
        self.encoder = nn.Sequential(
        nn.Linear(28 * 28, 64),
        nn.ReLU(),
        nn.Linear(64, 3))
        self.decoder = nn.Sequential(
        nn.Linear(3, 64),
        nn.ReLU(),
        nn.Linear(64, 28 * 28))

    def forward(self, x):
        embedding = self.encoder(x)
        return embedding

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=1e-3)
        return optimizer

    def training_step(self, train_batch, batch_idx):
        x, y = train_batch
        x = x.view(x.size(0), -1)
        z = self.encoder(x)    
        x_hat = self.decoder(z)
        loss = F.mse_loss(x_hat, x)
        self.log('train_loss', loss)
        return loss

    def validation_step(self, val_batch, batch_idx):
        x, y = val_batch
        x = x.view(x.size(0), -1)
        z = self.encoder(x)
        x_hat = self.decoder(z)
        loss = F.mse_loss(x_hat, x)
        self.log('val_loss', loss)

# data
dataset = MNIST('data', train=True, download=False, transform=transforms.ToTensor())
mnist_train, mnist_val = random_split(dataset, [55000, 5000])

train_kwargs = {'shuffle': True}
val_kwargs = {'shuffle': False}
cuda_kwargs = {'num_workers': int(os.environ["SLURM_CPUS_PER_TASK"]), 'pin_memory': True}
train_kwargs.update(cuda_kwargs)
val_kwargs.update(cuda_kwargs)

train_loader = DataLoader(mnist_train, batch_size=32, **train_kwargs)
val_loader = DataLoader(mnist_val, batch_size=32, **val_kwargs)

# model
model = LitAutoEncoder()

# training
trainer = pl.Trainer(gpus=1, num_nodes=1, precision=32, limit_train_batches=0.5, enable_progress_bar=False, max_epochs=10)
trainer.fit(model, train_loader, val_loader)
```

## Multi-GPU Example

Let's work through this [example](https://pytorch-lightning.readthedocs.io/en/latest/notebooks/lightning_examples/cifar10-baseline.html) where a modified ResNet-18 model is trained on [CIFAR-10](https://en.wikipedia.org/wiki/CIFAR-10). Here is the application script:

```python
import os

import torch
import torch.nn as nn
import torch.nn.functional as F
import torchvision
from pl_bolts.datamodules import CIFAR10DataModule
from pl_bolts.transforms.dataset_normalizations import cifar10_normalization
from pytorch_lightning import LightningModule, Trainer, seed_everything
from pytorch_lightning.callbacks import LearningRateMonitor
from pytorch_lightning.loggers import TensorBoardLogger
from torch.optim.lr_scheduler import OneCycleLR
from torch.optim.swa_utils import AveragedModel, update_bn
from torchmetrics.functional import accuracy

seed_everything(7)

BATCH_SIZE = 256
PATH_DATASETS = os.environ.get("PATH_DATASETS", ".")
NUM_WORKERS = int(os.environ["SLURM_CPUS_PER_TASK"])
NUM_NODES = int(os.environ["SLURM_NNODES"])
ALLOCATED_GPUS_PER_NODE = int(os.environ["SLURM_GPUS_ON_NODE"])

train_transforms = torchvision.transforms.Compose(
    [
        torchvision.transforms.RandomCrop(32, padding=4),
        torchvision.transforms.RandomHorizontalFlip(),
        torchvision.transforms.ToTensor(),
        cifar10_normalization(),
    ]
)

test_transforms = torchvision.transforms.Compose(
    [
        torchvision.transforms.ToTensor(),
        cifar10_normalization(),
    ]
)

cifar10_dm = CIFAR10DataModule(
    data_dir=PATH_DATASETS,
    batch_size=BATCH_SIZE,
    num_workers=NUM_WORKERS,
    train_transforms=train_transforms,
    test_transforms=test_transforms,
    val_transforms=test_transforms,
)

def create_model():
    model = torchvision.models.resnet18(pretrained=False, num_classes=10)
    model.conv1 = nn.Conv2d(3, 64, kernel_size=(3, 3), stride=(1, 1), padding=(1, 1), bias=False)
    model.maxpool = nn.Identity()
    return model

class LitResnet(LightningModule):
    def __init__(self, lr=0.05):
        super().__init__()

        self.save_hyperparameters()
        self.model = create_model()

    def forward(self, x):
        out = self.model(x)
        return F.log_softmax(out, dim=1)

    def training_step(self, batch, batch_idx):
        x, y = batch
        logits = self(x)
        loss = F.nll_loss(logits, y)
        self.log("train_loss", loss)
        return loss

    def evaluate(self, batch, stage=None):
        x, y = batch
        logits = self(x)
        loss = F.nll_loss(logits, y)
        preds = torch.argmax(logits, dim=1)
        acc = accuracy(preds, y)

        if stage:
            self.log(f"{stage}_loss", loss, prog_bar=True)
            self.log(f"{stage}_acc", acc, prog_bar=True)

    def validation_step(self, batch, batch_idx):
        self.evaluate(batch, "val")

    def test_step(self, batch, batch_idx):
        self.evaluate(batch, "test")

    def configure_optimizers(self):
        optimizer = torch.optim.SGD(
            self.parameters(),
            lr=self.hparams.lr,
            momentum=0.9,
            weight_decay=5e-4,
        )
        steps_per_epoch = 45000 // BATCH_SIZE
        scheduler_dict = {
            "scheduler": OneCycleLR(
                optimizer,
                0.1,
                epochs=self.trainer.max_epochs,
                steps_per_epoch=steps_per_epoch,
            ),
            "interval": "step",
        }
        return {"optimizer": optimizer, "lr_scheduler": scheduler_dict}

model = LitResnet(lr=0.05)
model.datamodule = cifar10_dm

trainer = Trainer(
    gpus=ALLOCATED_GPUS_PER_NODE,
    num_nodes=NUM_NODES,
    strategy='ddp',
    precision=32,
    max_epochs=10,
    progress_bar_refresh_rate=10,
    logger=TensorBoardLogger("lightning_logs/", name="resnet"),
    callbacks=[LearningRateMonitor(logging_interval="step")],
)

trainer.fit(model, cifar10_dm)
trainer.test(model, datamodule=cifar10_dm)
```


### Step 1: Download the data

The compute nodes do not have internet access so download the data on the login node:

```
$ cd multi_gpu_training/03_pytorch_lightning/multi
(bolts) $ python download_cifar10.py
```

### Step 2: Submit the Job

Below is the Slurm script:

([see script](job.slurm)).

If you want to change the Pytorch distributed backend you can 

```bash
export PL_TORCH_DISTRIBUTED_BACKEND=gloo
```

Submit the job:

```bash
$ sbatch job.slurm  # edit your email address in job.slurm before submitting
```

How does the training time decrease in going from 1 to 2 to 4 GPUs? What happens if you use `precision=16`?

## Numerical Precision

You can try adjusting the `precision` to accelerate training. 

NOTE: The choice of `precision="bf16"` can only be used with PyTorch 1.10.

## Debugging

For troubleshooting NCCL try adding these environment variables to your Slurm script:

```
export NCCL_DEBUG=INFO
export NCCL_DEBUG_SUBSYS=ALL
```

## Useful Links

+ [PyTorch Lightning and Slurm](https://pytorch-lightning.readthedocs.io/en/latest/clouds/cluster.html)  
+ [PyTorch LightningLite](https://pytorch-lightning.readthedocs.io/en/latest/starter/lightning_lite.html)
