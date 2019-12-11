import numpy as np
import math as m

def get_distance(coordinate1, coordinate2):
    a = coordinate1[0] - coordinate2[0]
    b = coordinate1[1] - coordinate2[1]
    return m.sqrt(a**2 + b**2)

def distances(lst):
    a = []
    for i in range(len(lst)):
        a.append(get_distance(lst[i], lst[(i+1)%len(lst)]))
    return a

def get_Chord(x,root,taper,span):
    return (root - ((x/span)*(1-taper)*root))

def PolyArea(x,y):
        return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

class halfwing():
    def __init__(self, span, taper, rootChord, thickness, nstringer, Astringer):
        self.span = span
        self.wingBoxArray = np.array([[0.15,-0.0506],
                            [ 0.6,   -0.0444],
                            [ 0.6,     0.0541],
                            [ 0.15,    0.0505]])
        self.taper = taper
        self.rootChord = rootChord
        self.tipChord = self.rootChord * self.taper
        self.thickness = thickness
        self.Astringer = Astringer
        self.nstringer = nstringer
        self.perimeter = None
        self.iterWingBox = None
        self.chord = None

    @staticmethod
    def get_Chord(x,root,taper,span):
        return (root - ((x/span)*(1-taper)*root))

    def getChord(self, x):
        self.chord = get_Chord(x, self.rootChord,self.taper,self.span)

    def getPerimeter(self, x):
        self.getChord(x)
        self.iterWingBox = self.wingBoxArray
        self.iterWingBox = self.iterWingBox * self.chord
        self.perimeter = sum(distances(list(self.iterWingBox)))

    def getVolume(self):
        area = (self.getTotalArea() * self.thickness) + self.getStringerArea(self.nstringer, self.Astringer)
        return area

    def getStringerArea(self, nstringer, Astringer):
        stringerarea = self.span * nstringer * Astringer
        return stringerarea

    def getTotalArea(self):
        dx = 0.001
        x = 0
        totalarea = 0
        while x < self.span:
            self.getPerimeter(x)
            totalarea += dx * self.perimeter
            x += dx
        return totalarea
    
    def getAreaAtSpan(self, span):
        self.getPerimeter(span)
        return self.thickness * self.perimeter

stringerHeight = 0.05
stringerWidth = 0.02
stringerThickness = 0.01
numberofstringers = 30
wingboxthickness = 0.002

def stringerArea(height, width, thickness):
    return (height + width) * thickness
areaofstringer = stringerArea(stringerHeight,stringerWidth,stringerThickness)

ourPlane = halfwing(21.74/2,0.4,3.45, wingboxthickness, numberofstringers, areaofstringer)


print(ourPlane.getAreaAtSpan(4)) #Change this value in m
