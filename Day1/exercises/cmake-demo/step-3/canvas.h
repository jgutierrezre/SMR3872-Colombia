// -*- c++ -*-

#ifndef CANVAS_H
#define CANVAS_H

#include "color.h"

class Canvas
{
public:
    Canvas(int width, int height, const Color &background);
    virtual ~Canvas();

public:
    void clear(const Color &background);
    void fill(int x1, int y1, int x2, int y2, const Color &forground);

    void write_ppm(const char *filename);
    void write_jpg(const char *filename);

  private:
    int width, height;
    Color *buffer;

  private:
    // disable default and copy constructor
    Canvas() {};
    Canvas(const Canvas &) {};
};


#endif
