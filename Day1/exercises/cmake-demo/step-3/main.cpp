// -*- c++ -*-

#include "canvas.h"

int main(int, char**)
{
    Canvas *flag = new Canvas(640,480,white);
    flag->write_ppm("white.ppm");

    flag->fill(0,0,213,480,red);
    flag->fill(426,0,640,480,green);
    flag->write_ppm("it_flag.ppm");
    flag->write_jpg("it_flag.jpg");

    flag->clear(black);
    flag->fill(0,320,640,480,yellow);
    flag->fill(0,160,640,320,red);
    flag->write_ppm("de_flag.ppm");
    flag->write_jpg("de_flag.jpg");

    delete flag;
    return 0;
}
