import rhinoscriptsyntax as rs

fileID1 = raw_input("File Name")
count1 = 1
patches = rs.GetObjects("Select Surfaces to Export",rs.filter.surface)
for patch in patches:
    fileID = fileID1+str(count1)
    fileText = open(fileID,"w")
    count2 = 0
    knotCount = rs.SurfaceKnotCount(patch)
    knots = rs.SurfaceKnots(patch)
    pointCount = rs.SurfacePointCount(patch)
    points = rs.SurfacePoints(patch)
    weights = rs.SurfaceWeights(patch)
    for knot in knots:
        fileText.write(str(knotCount[count2]+2))
        fileText.write("\n")
        difference = 0.0
        fileText.write(str(knot[0]))
        fileText.write(" ")
        for i in range(0,knotCount[count2]):
            input = knot[i]
            fileText.write("%.6f" % (input))
            fileText.write(" ")
        fileText.write(str(knot[-1]))
        fileText.write("\n")
        count2 += 1
    fileText.write(str(pointCount[0]*pointCount[1]))
    fileText.write("\n")
    i=0
    for v in range(pointCount[0]):
        for u in range(pointCount[1]):
            fileText.write("%.6f" % (points[i][0]))
            fileText.write(" ")
            fileText.write("%.6f" % (points[i][1]))
            fileText.write(" ")
            fileText.write("%.6f" % (points[i][2]))
            fileText.write(" ")
            fileText.write("%.6f" % (weights[i]))
            fileText.write("\n")
            i += 1
    count1 += 1
    fileText.close()
