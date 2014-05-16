'''
Created on Oct 8, 2010

@author: quilan
'''

#========================================

class NewtonsMethod(object):
    def __init__(self, coefficientTerms, powerTerms):
        self.ACC=10.0**10   #When we consider having reached a zero. if |f(z)| < 1/ACC, it's a zero.
        self.LOOPMAX=100    #How far we calculate before bailing.
        self.initTerms(coefficientTerms, powerTerms)

    def initTerms(self, coefficientTerms, powerTerms):
        self.powerTerms, self.coefficientTerms = self.consolidateTerms(coefficientTerms, powerTerms)
        self.formula = self.getFormula(self.powerTerms, self.coefficientTerms)
        self.degree = max(self.powerTerms)
        self.f = lambda z: sum(c*z**p for c,p in zip(self.coefficientTerms,self.powerTerms))
        self.df = lambda z: sum(c*p*z**(p-1) for c,p in zip(self.coefficientTerms,self.powerTerms))

    #Takes a series of random powers/coefficients & combines like powers.
    #Basically, doing: a*z^n + b*z^n = (a+b)*z^n
    def consolidateTerms(self,cterm,pterm):
        eq={}
        for a,b in zip(cterm,pterm):
            try: eq[b]+=a
            except: eq[b]=a

        for x in eq.keys():
            if eq[x]==0: del eq[x]

        if len(eq.items()) == 0: eq = { 0:0 }
        pterm,cterm=zip(*sorted(eq.items(),reverse=True))
        return pterm,cterm

    #Creates a human-readable representation of the powers/coefficients
    def getFormula(self,pterm,cterm):
        formula=""
        for p,c in zip(pterm,cterm):
            if p==1: pst="z"
            elif p==0: pst=""
            else: pst="z^%s"%p

            if c>1: cst="+%s"%c
            elif c==1 and p>0: cst="+"
            elif c==1 and p==0: cst="+1"
            elif c==-1 and p>0: cst="-"
            else: cst="%s"%c

            formula="%s%s%s"%(formula,cst,pst)

        if formula[0]=='+': formula=formula[1:]
        return formula

    #Perform Newton's method from the starting point z
    def newtonIter(self,z):
        try:
            for iloop in xrange(self.LOOPMAX):
                z -= self.f(z)/self.df(z)
                if self.ACC*abs(self.f(z))<1: break
            else: iloop=self.LOOPMAX
        except: iloop=self.LOOPMAX

        return z,iloop

    #Performs the Durand-Kerner method of finding all roots of a function
    #See the wiki for an idea on how this method works
    def getAllRoots(self):

        #Let's go for a solid 10 attempts before giving up
        for attempt in xrange(10):
            try:
                #We want to find a random complex number that is not a root of unity
                while True:
                    root = complex(random(), random())
                    if self.f(root) == 0: continue
                    break

                #Initialize each initial root to a fairly random number
                roots = [root**i for i in xrange(self.degree)]

                #Keep iterating until we find 
                for _ in xrange(self.LOOPMAX):
                    newRoots = roots[:]
                    for i in xrange(self.degree):
                        z = roots[i]
                        denominator = reduce(operator.mul, ( z-roots[j] for j in xrange(self.degree) if i!=j ), 1 );
                        newRoots[i] = z - self.f(z)/denominator
                    roots = newRoots

                    #If every single root is less than 1/ACC, bail because we've found them all!
                    if all( self.ACC*abs(self.f(z)) < 1 for z in roots ):
                        return roots

            #If we had some weird divide-by-zero exception or something, restart with a new root seed
            except: continue

        #If we've failed all attempts at finding the roots, we're boned, so just bail
        #We'll pick the roots up later in the drawing process
        return [];

class NewtonDrawer(object):
    def __init__(self,equation,imageWidth,imageHeight,rect):
        self.equation = equation
        self.imageWidth = imageWidth
        self.imageHeight = imageHeight
        self.rect = rect[:]
        self.data = [ [ () for _ in xrange(imageWidth) ] for _ in xrange(imageHeight) ]
        self.STDDEV=20.0    #Shading. Higher == bright+low contrast, Lower == dark, more contrast

    #Returns (r,g,b) from a hue from 0 to 360.
    #0 degrees == blue
    #120 degrees == red
    #240 degrees == green
    def getColor(self,hue):
        if hue<0: hue+=360
        if hue<120: g=0; r=hue; b=120-hue;
        elif hue<240: b=0; g=hue-120; r=240-hue;
        else: r=0; b=hue-240; g=360-hue;

        r/=120.0; g/=120.0; b/=120.0;
        m=max(r,g,b)/255.0
        r/=m; g/=m; b/=m;

        return [int(r),int(g),int(b)]

    #Saves data to a binary file
    def saveDataToFile(self, fileName):
        with open(fileName, "wb") as f:
            f.write(struct.pack("i", len(self.equation.coefficientTerms)))
            for c,p in zip(self.equation.coefficientTerms, self.equation.powerTerms):
                f.write(struct.pack("ii", c, p))

            f.write(struct.pack("II", self.imageWidth, self.imageHeight))
            f.write(struct.pack("dddd", *self.rect))

            for y in xrange(self.imageWidth):
                for x in xrange(self.imageHeight):
                    z,iloop = self.data[y][x]
                    f.write(struct.pack("ddI", z.real, z.imag, iloop))

    #Reads data from a binary file
    def getDataFromFile(self, fileName):
        with open(fileName, "rb") as f:
            cterms = []; pterms = []
            s = "i"; sz = struct.calcsize(s)
            count, = struct.unpack(s, f.read(sz))
            for _ in xrange(count):
                s = "ii"; sz = struct.calcsize(s)
                c,p = struct.unpack(s, f.read(sz))
                cterms.append(c)
                pterms.append(p)
            self.equation.initTerms(cterms,pterms)

            s = "II"; sz = struct.calcsize(s)
            self.imageWidth, self.imageHeight = struct.unpack(s, f.read(sz))
            self.data = [ [ () for _ in xrange(self.imageWidth) ] for _ in xrange(self.imageHeight) ]

            s = "dddd"; sz = struct.calcsize(s)
            self.rect = struct.unpack(s, f.read(sz))

            s = "ddI"; sz = struct.calcsize(s)
            for y in xrange(self.imageWidth):
                for x in xrange(self.imageHeight):
                    zr,zi,iloop = struct.unpack(s, f.read(sz))
                    self.data[y][x] = complex(zr,zi), iloop

    #Maps an x,y coordinate to a complex number
    def pixelToComplex(self):
        realOrigin = self.rect[0]
        imagOrigin = self.rect[1]
        realRange  = self.rect[2]
        imagRange  = self.rect[3]

        def scale(x,y):
            sx = realOrigin + 1.0 * realRange * x / self.imageWidth
            sy = imagOrigin + 1.0 * imagRange * y / self.imageHeight
            return sx,sy
        return scale

    #Maps a complex number to an x,y coordinate
    def complexToPixel(self):
        realOrigin = self.rect[0]
        imagOrigin = self.rect[1]
        realRange  = self.rect[2]
        imagRange  = self.rect[3]

        def unscale(sx,sy):
            x = 1.0 * (sx-realOrigin) * self.imageWidth  / realRange
            y = 1.0 * (sy-imagOrigin) * self.imageHeight / imagRange
            return x,y
        return unscale

    #Draw the fractal on the images
    def renderFractal(self):
        pixelToComplex = self.pixelToComplex()

        lst=time()
        for y in xrange(self.imageHeight):
            if time()-lst>300: lst=time(); print "\t%d%% finished rendering"%(y*100/self.imageHeight)
            for x in xrange(self.imageWidth):
                #Create the initial z-value, and then perform Newton's Method on it.
                r,i=pixelToComplex(x,y)
                self.data[y][x]=self.equation.newtonIter(complex(r,i))

    #Draw the fractal onto an image
    def drawFractal(self, drawSurface):
        import Image, ImageDraw
        roots = { r:self.getColor(atan2(r.imag,r.real)*180/pi) for r in self.equation.getAllRoots() }

        lst=time()
        for y in xrange(self.imageHeight):
            if time()-lst>300: lst=time(); print "\t%d%% finished rendering"%(y*100/self.imageHeight)
            for x in xrange(self.imageWidth):
                z,iloop = self.data[y][x]

                s=0
                c=[0,0,0]
                if(iloop<self.equation.LOOPMAX):
                    mn,mr=0.05,z
                    #Find the nearest existing root within 0.05 of our final z value.
                    for r,c in roots.items():
                        if(abs(r-z)<mn):
                            mn=abs(r-z)
                            mr=r

                    #If we couldn't find an existing root within 0.05 of our z value, make a new one!
                    if(mr not in roots):
                        #print "\tNew root: %s\t\t\t%s"%(mr,self.equation.f(z))
                        roots[mr]=self.getColor(atan2(mr.imag, mr.real)*180/pi)

                    c=roots[mr]

                    ##Should be a value between 0 and 1, used to smooth between loop counts
                    ff=(self.equation.ACC*abs(self.equation.f(z)))**0.1
                    #Basic gaussian curve, e^(-x^2/o^2)
                    s=exp(-((iloop+ff)/self.STDDEV)**2)
                    s=int(s*255)
                drawSurface.point((x,y),tuple([int(v*s/255.0) for v in c]))

        complexToPixel = self.complexToPixel();
        #Draw a little white circle around each found root
        for rt in roots.keys():
            r,i=rt.real,rt.imag
            x,y = complexToPixel(r,i)
            drawSurface.ellipse((x-2,y-2,x+2,y+2),fill=(255,255,255))

#========================================

def drawNewtonFractals():
    imageWidth=imageHeight=200
    imageRect=[-1.5,-1.5,3,3]

    #Previously found formulas that are rendered first, before doing random ones
    saved = []
    saved.append([[13,6,1,0],[1,-3,1,-1]])      #1*z^13 - 3*z^6 + 1*z - 1
    saved.append([[12,4,1,0],[2,-2,-2,-1]])     #2*z^12 - 2*z^4 - 2*z - 1
    saved.append([[7,1,0],[2,-1,2]])            #2*z^7 - 1*z + 2
    saved.append([[5,3,1,0],[1,3,1,3]])         #1*z^5 + 3*z^3 + 1*z + 3
    saved.append([[6,4,2,0],[1,2,-3,1]])        #etc.
    saved.append([[4,2,0],[2,-2,-1]])
    saved.append([[6,2,0],[2,2,2]])
    saved.append([[7,4,1],[2,3,2]])
    saved.append([[12,3,0],[2,-2,1]])
    saved.append([[3,2],[2,-2]])

    #When creating random formulas, each term will be (-MC .. MC) * z^(0..MP)
    MP=15
    MC=4

    cache=set()
    for iteration in xrange(10000):
        #If we're still loading saved formulas, use the saved variables
        if iteration<len(saved):
            pterm,cterm = saved[iteration]
        else:
            #else, create a random formula to render
            nterm=randint(2,4)
            cterm = [randint(-MC,MC) for _ in xrange(nterm)]
            pterm = [randint(0,MP) for _ in xrange(nterm)]
            pterm.append(0);
            cterm.append(randint(-MC,MC));

        equation = NewtonsMethod(cterm, pterm)
        newtonDrawer = NewtonDrawer(equation,imageWidth,imageHeight,imageRect)

        if equation.formula in cache: continue
        cache.add(equation.formula)

        print "Formula: %s"%equation.formula

        start=time()
        newtonDrawer.renderFractal()
        print "Total render time: %fs"%(time()-start)

        start=time()
        imgRoots = Image.new("RGB", (imageWidth, imageHeight), "#FFFFFF")
        newtonDrawer.drawFractal(ImageDraw.Draw(imgRoots))
        print "Total draw time: %fs"%(time()-start)

        imgRoots.save("pics\%s.png"%equation.formula)
        print

#========================================

from math import exp, atan2, pi
from random import randint, random
from time import time
from os import mkdir
import operator
import struct

if __name__ == "__main__":
    try: mkdir("pics")
    except: pass
    drawNewtonFractals()

#========================================