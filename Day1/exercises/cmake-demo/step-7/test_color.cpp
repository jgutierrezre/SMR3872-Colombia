// -*- c++ -*-

#include "color.h"
#include "gtest/gtest.h"

TEST(Color, constructor)
{
        Color def;
        ASSERT_EQ(def.r,0);
        ASSERT_EQ(def.g,0);
        ASSERT_EQ(def.b,0);
        def.r = 10;
        def.g = 'c';
        Color cpy(def);
        ASSERT_EQ(cpy.r,10);
        ASSERT_EQ(cpy.g,'c');
        ASSERT_EQ(cpy.b,0);
}

TEST(Color, presets)
{
        ASSERT_EQ(white.r,255);
        ASSERT_EQ(white.g,255);
        ASSERT_EQ(white.b,255);
        ASSERT_EQ(black.r,0);
        ASSERT_EQ(black.g,0);
        ASSERT_EQ(black.b,0);
        ASSERT_EQ(red.r,255);
        ASSERT_EQ(red.g,0);
        ASSERT_EQ(red.b,0);
        ASSERT_EQ(green.r,0);
        ASSERT_EQ(green.g,255);
        ASSERT_EQ(green.b,0);
        ASSERT_EQ(blue.r,0);
        ASSERT_EQ(blue.g,0);
        ASSERT_EQ(blue.b,255);
        ASSERT_EQ(yellow.r,255);
        ASSERT_EQ(yellow.g,255);
        ASSERT_EQ(yellow.b,0);
}

