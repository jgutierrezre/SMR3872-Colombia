// -*- c++ -*-

#include "canvas.h"
#include <cstdio>

#if defined(HAVE_JPEG_LIB)
#include <jpeglib.h>
#endif

Canvas::Canvas(int _width, int _height, const Color &_bg)
        : width(_width), height(_height), buffer(nullptr)
{
    const int npixel = width*height;
    buffer = new Color[npixel];
    for (int i=0; i < npixel; ++i) buffer[i] = _bg;
}

Canvas::~Canvas()
{
    delete[] buffer;
}


void Canvas::clear(const Color &_bg)
{
    const int npixel = width*height;
    for (int i=0; i < npixel; ++i) buffer[i] = _bg;
}

void Canvas::fill(int x1, int y1, int x2, int y2, const Color &foreground)
{
    if (x1 < 0) x1 = 0;
    if (y1 < 0) y1 = 0;
    if (x2 >= width)  x2 = width-1;
    if (y2 >= height) y2 = height-1;

    for (int j=y1; j <= y2; ++j)
        for (int i=x1; i <= x2; ++i)
            buffer[j*width+i] = foreground;
}

void Canvas::write_ppm(const char *filename)
{
    FILE *fp = fopen(filename,"wb");
    if (!fp) return;

    fprintf(fp,"P6\n%d %d\n255\n",width,height);
    for (int i=0; i < height; ++i) {
        fwrite(&buffer[i*width].r,3,width,fp);
    }
    fclose(fp);
}

void Canvas::write_jpg(const char *filename)
{
#if defined(HAVE_JPEG_LIB)
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    JSAMPROW row_pointer;

    FILE *fp = fopen(filename,"wb");
    if (!fp) return;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo,fp);
    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo,85,TRUE);
    jpeg_start_compress(&cinfo,TRUE);

    while (cinfo.next_scanline < cinfo.image_height) {
        row_pointer = (JSAMPROW) &buffer[width*cinfo.next_scanline].r;
        jpeg_write_scanlines(&cinfo,&row_pointer,1);
    }
    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);
    fclose(fp);
#else
    fprintf(stderr,"No JPEG support available.\n");
#endif
}
