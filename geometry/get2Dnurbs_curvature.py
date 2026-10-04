import rhinoscriptsyntax as rs

fileID1 = raw_input("File Name")
count1 = 1
patches = rs.GetObjects("Select Surfaces to Export", rs.filter.surface)

for patch in patches:
    fileID = fileID1 + str(count1) + "_curvature.txt"
    fileText = open(fileID, "w")

    # Define the 5x5 sampling grid
    samples = 5
    domainU = rs.SurfaceDomain(patch, 0)
    domainV = rs.SurfaceDomain(patch, 1)

    fileText.write("U_Param, V_Param, K1, K2\n")

    for i in range(samples):
        for j in range(samples):
            # Calculate equally distributed parameters
            u = domainU[0] + (i / float(samples - 1)) * (domainU[1] - domainU[0])
            v = domainV[0] + (j / float(samples - 1)) * (domainV[1] - domainV[0])

            # Read curvature from Rhino
            curvatureData = rs.SurfaceCurvature(patch, [u, v])

            if curvatureData:
                k1 = curvatureData[2] # Maximum principal curvature
                k2 = curvatureData[4] # Minimum principal curvature
                fileText.write("%.4f" % (u))
                fileText.write(" ")
                fileText.write("%.4f" % (v))
                fileText.write(" ")
                fileText.write("%.4f" % (k1))
                fileText.write(" ")
                fileText.write("%.4f" % (k2))

    fileText.close()
    count1 += 1
