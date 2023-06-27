import numpy as np
import matplotlib.pyplot as plt

#c1 = (255.0, 255.0, 255.0)
c2 = (255.0, 0.0, 0.0)
c3 = (0.0, 255.0, 0.0)
c4 = (0.0, 0.0, 255.0)

def get_color(col_tup):
    return "#{0:02X}{1:02X}{2:02X}".format(int(col_tup[0]), int(col_tup[1]), int(col_tup[2]))

print()
print("const uint8_t color_lut[256][3] = {")

img = np.zeros((1,256,3))

# for i in range(128):
#     r = i/85 * (c2[0]-c1[0]) + c1[0]
#     g = i/85 * (c2[1]-c1[1]) + c1[1]
#     b = i/85 * (c2[2]-c1[2]) + c1[2]
#     img[0][i][0] = int(r)
#     img[0][i][1] = int(g)
#     img[0][i][2] = int(b)
#     print("\t{"+str(int(r))+","+str(int(g))+","+str(int(b))+"},")

for i in range(128):
    r = i/128 * (c3[0]-c2[0]) + c2[0]
    g = i/128 * (c3[1]-c2[1]) + c2[1]
    b = i/128 * (c3[2]-c2[2]) + c2[2]
    img[0][i][0] = int(r)
    img[0][i][1] = int(g)
    img[0][i][2] = int(b)
    print("\t{"+str(int(r))+","+str(int(g))+","+str(int(b))+"},")

for i in range(128):
    r = i/128 * (c4[0]-c3[0]) + c3[0]
    g = i/128 * (c4[1]-c3[1]) + c3[1]
    b = i/128 * (c4[2]-c3[2]) + c3[2]
    img[0][i+128][0] = int(r)
    img[0][i+128][1] = int(g)
    img[0][i+128][2] = int(b)
    print("\t{"+str(int(r))+","+str(int(g))+","+str(int(b))+"},")

print("}")
print()

img = img.astype(int)
plt.imshow(img, aspect='auto')
plt.show()