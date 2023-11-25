import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
tree = ET.parse('arch.geo')
 
root = tree.getroot()
jointlist=root[0][1]
n=jointlist.__len__()
for i in range(n):
    value=jointlist[i].text
    values=value.split("     ")
    print(values)
    plt.plot([float(values[0]),float(values[2])],[float(values[1]),float(values[3])])
plt.show()

