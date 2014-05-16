NewtonsMethod
=============
tl-dr: This program will create a 'pics' directory where the program is
called, and fill it with fractals created by Newton's Method.

Here's an old python script I wrote to generate random Newton's Method fractal
patterns. It's hobbled together with some duct tape & string, so don't expect
much in terms of modern pretty code.

To play around with things, edit the "drawNewtonFractals" function. The key
things to potentially change are the "imageWidth"/"imageHeight" variables, the
"imageRect" coordinates (minimum complex point, real width, imag height) and
finally the "saved" array. These are just some fun looking formulas I noticed
when combing through hundreds of pictures after a day or so of running. If you
want to do nothing but spawn random images, comment out the entire saved array
(minus the declaration, of course).

===Requirements===

This does require the Python Imaging Library:
http://www.pythonware.com/products/pil/

If you can get PIL to work with PyPy, you'll be in rendering heaven as it
should be nice & zippy. If not, things might take a bit; take a break, have
a nice cup of tea and enjoy the world.

=============