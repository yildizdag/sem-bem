import rhinoscriptsyntax as rs

fileID1 = raw_input("File Name")
count1 = 1
patches = rs.GetObjects("Select Curves to Export",rs.filter.curve)
for patch in patches:
    fileID = fileID1+str(count1)
    fileText = open(fileID,"w")

    knotCount = rs.CurveKnotCount(patch)
    knot = rs.CurveKnots(patch)
    pointCount = rs.CurvePointCount(patch)
    points = rs.CurvePoints(patch)
    weights = rs.CurveWeights(patch)
    fileText.write(str(knotCount+2))
    fileText.write("\n")
    norm = knot[-1]
    fileText.write(str(0))
    fileText.write(" ")
    for i in range(0,knotCount):
        fileText.write("%.8f" % (knot[i]/norm))
        fileText.write(" ")
    fileText.write(str(1))
    fileText.write("\n")

    fileText.write(str(pointCount))
    fileText.write("\n")
    i=0
    for v in range(pointCount):
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
