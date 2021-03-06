#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division, print_function
import svgwrite

__doc__ = """
The internal plotting library. The main purpose of this library
is to give a graphical representation of the dna chromatogram
data. Still this library may also be util/usefull for debugging
purpose.
"""

# the different colors
colorCodes = [
    (255,0,0),
    (0,255,0),
    (0,0,255),
    (0,0,0),
    (0,255,255),
    (255,0,255),
    (255,255,0),
    (200,200,200)
]

def _unicode(obj):
    """
        Needed so that code works with both python2 and pyhton3.
    """
    try: return unicode(obj)
    except: return str(obj)

def plotChromatogram(filename, chromatogram, keys=["A","C","T","G"]):
    """
        Plot chromatogram data to a file.

        @param filename      The output file (*.svg) -
        the file the plot should get saved to.
        @param chromatogram  The chromatogram data to plot.
        @param keys          The keys of the traces in the chromatogram to plot.
    """
    # settings - maybe need to adjust these values later ...
    offsetX, offsetY = 50, 50
    # get the maximal value of the chromatogram
    maxChromVal = float("-inf")
    for key in keys:
        maxChromVal = max(maxChromVal, max(chromatogram[key]))
    # get the length of the chromatogram
    chromLength = 0
    for key in keys:
        chromLength = max(chromLength, len(chromatogram[key]))
    # create an svg object
    svg = svgwrite.Drawing(filename=filename,
        size=(2 * offsetX + chromLength, 2 * offsetY + maxChromVal))
    # plot surrounding box
    drawbox = svg.rect(insert=(offsetX, offsetY), size=(chromLength, maxChromVal),
        fill="white", stroke="black")
    svg.add(drawbox)
    # plot legend
    for i, key in enumerate(keys):
        # add the box
        style = "fill:rgba(%s,%s,%s,0.3);" % colorCodes[i % len(colorCodes)]
        style += "stroke:rgb(%s,%s,%s)" % colorCodes[i % len(colorCodes)]
        box = svg.rect(insert=(offsetX + 15, offsetY + 15 + i * 24), size=(24, 12), style=style)
        svg.add(box)
        # add the text
        style = "font-size:20px"
        text = svg.text(_unicode(key), insert=(offsetX + 52, offsetY + 27 + i * 24), style=style)
        svg.add(text)
    # plot bottom scale/ruler
    for i in range(0, chromLength, 100):
        # ruler line
        line = svg.line(
            start=(offsetX + i, offsetY + maxChromVal),
            end=(offsetX + i, offsetY + maxChromVal + 10),
            style="stroke:black;stroke-width:1"
        )
        svg.add(line)
        # ruler text
        text = svg.text(str(i), insert=(offsetX + i - 5, offsetY + maxChromVal + 28))
        svg.add(text)
    # plot each trace
    for i, key in enumerate(keys):
        # plot the trace
        points = [(offsetX, offsetY + maxChromVal)] # starting point is always "0", "0"
        lastPoint = offsetX
        for j, value in enumerate(chromatogram[key]):
            lastPoint = offsetX + j
            points.append((lastPoint, offsetY + maxChromVal - value))
        points.append((lastPoint, offsetY + maxChromVal)) # ending point is always "end", "0"
        # create the polygon
        style = "fill:rgba(%s,%s,%s,0.3);" % colorCodes[i % len(colorCodes)]
        style += "stroke:rgb(%s,%s,%s)" % colorCodes[i % len(colorCodes)]
        poly = svg.polygon(points, style=style)
        # add the polygon to the figure
        svg.add(poly)
    # save and return the svg object
    svg.save()

def plotCwt(filename, cwt, keys=["A","C","T","G"]):
    """
        Plot the continuous wavelet transformation data to a file.

        @param filename      The output file (*.svg) -
        the file the plot should get saved to.
        @param cwt           The continuous wavelet transformation data.
        @param keys          The keys of the traces in the cwt to plot.
    """
    # settings - maybe need to adjust these values later ...
    offsetX, offsetY = 50, 50
    # get the maximal and minimal value of the cwt (as well as the length)
    cwtLength, maxVal, minVal = 0, float("-inf"), float("inf")
    for key in keys:
        maxVal = max(maxVal, max(cwt[key]))
        minVal = min(minVal, min(cwt[key]))
        cwtLength = max(cwtLength, len(cwt[key]))
    # create an svg object
    svg = svgwrite.Drawing(filename=filename,
        size=(2 * offsetX + cwtLength, 2 * offsetY + maxVal - minVal))
    # plot surrounding box
    drawbox = svg.rect(insert=(offsetX, offsetY), size=(cwtLength, maxVal - minVal),
        fill="white", stroke="black")
    svg.add(drawbox)
    # plot legend
    for i, key in enumerate(keys):
        # add the box
        style = "fill:rgba(%s,%s,%s,0.3);" % colorCodes[i % len(colorCodes)]
        style += "stroke:rgb(%s,%s,%s)" % colorCodes[i % len(colorCodes)]
        box = svg.rect(insert=(offsetX + 15, offsetY + 15 + i * 24), size=(24, 12), style=style)
        svg.add(box)
        # add the text
        style = "font-size:20px"
        text = svg.text(_unicode(key), insert=(offsetX + 52, offsetY + 27 + i * 24), style=style)
        svg.add(text)
    # plot bottom scale/ruler
    for i in range(0, cwtLength, 100):
        # ruler line
        line = svg.line(
            start=(offsetX + i, offsetY + maxVal - minVal),
            end=(offsetX + i, offsetY + maxVal - minVal + 10),
            style="stroke:black;stroke-width:1"
        )
        svg.add(line)
        # ruler text
        text = svg.text(str(i), insert=(offsetX + i - 5, offsetY + maxVal - minVal + 28))
        svg.add(text)
    # plot each trace
    for i, key in enumerate(keys):
        # plot the trace
        points = []
        for j, value in enumerate(cwt[key]):
            points.append((offsetX + j, offsetY + maxVal - value))
        # create the polygon
        style = "fill:rgba(0,0,0,0);"
        style += "stroke:rgb(%s,%s,%s)" % colorCodes[i % len(colorCodes)]
        poly = svg.polyline(points, style=style)
        # add the polygon to the figure
        svg.add(poly)
    # save and return the svg object
    svg.save()
