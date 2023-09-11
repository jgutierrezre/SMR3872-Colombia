#include <ctype.h>
#include <string.h>

char *makeupper(char *text)
{
    int len=strlen(text);
    char delta = 'A' - 'a';
    for (int i=0; i < len; ++i) {
        if (islower(text[i]))
            text[i] += delta;
    }
    return text;
}


