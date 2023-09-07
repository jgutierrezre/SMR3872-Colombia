// -*- c++ -*-

#ifndef COLOR_H
#define COLOR_H

struct Color
{
    Color() : r(0), g(0), b(0) {};
    Color(int _r, int _g, int _b)
            : r(_r), g(_g), b(_b) {};

    unsigned char r,g,b;
};

const Color black(0,0,0);
const Color white(255,255,255);
const Color red(255,0,0);
const Color green(0,255,0);
const Color blue(0,0,255);
const Color yellow(255,255,0);
const Color magenta(255,0,255);
const Color cyan(0,255,255);

#endif
