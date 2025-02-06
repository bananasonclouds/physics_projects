input_text = """
Date 1:
X Coordinate: 2.3612069587858016 AU
Y Coordinate: 1.6191920596530622 AU
Z Coordinate: -0.2764542590196006 AU
-------------------
Date 2:
X Coordinate: 2.6023272842677736 AU
Y Coordinate: 2.043247622471402 AU
Z Coordinate: -0.09263422488680328 AU
-------------------
Date 3:
X Coordinate: 2.5818084956596286 AU
Y Coordinate: 2.509502653300652 AU
Z Coordinate: 0.10948136513306542 AU
-------------------
Date 4:
X Coordinate: 2.297899431130209 AU
Y Coordinate: 2.9196309971394543 AU
Z Coordinate: 0.2872697246918256 AU
-------------------
Date 5:
X Coordinate: 1.8386780627328672 AU
Y Coordinate: 3.13647913923454 AU
Z Coordinate: 0.3812743272474866 AU
"""

# Split the input text into individual sections based on "Date"
sections = input_text.split("-------------------\n")

# Initialize a list to store the converted coordinates
converted_coordinates = []

# Iterate through the sections and extract numerical values
for section in sections:
    lines = section.strip().split('\n')
    x = float(lines[1].split(': ')[1].split(' AU')[0])
    y = float(lines[2].split(': ')[1].split(' AU')[0])
    z = float(lines[3].split(': ')[1].split(' AU')[0])
    converted_coordinates.append((x, y, z))

# Print the converted coordinates
for coordinate in converted_coordinates:
    print(coordinate)
