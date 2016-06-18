
import xml.etree.cElementTree as elementTree
import re
from matplotlib import pyplot as plt
import numpy as np
from svg.path import parse_path
from functools import lru_cache


# Flags for printing the output to the console.
PRINT_EXTRA_OUTPUT = False
PRINT_ELECTRODES = False

# Flags for logging output in data files.
LOG_EXTRA_OUTPUT = False
LOG_ELECTRODES = False

# Create an iterator over the trap image file
trap_XML_iterator = elementTree.iterparse('resources/RS1096.svg')
# Strip all namespaces
for _, element in trap_XML_iterator:
    if '}' in element.tag:
        element.tag = element.tag.split('}', 1)[1]
# Grab the root node of the iterator
root = trap_XML_iterator.root
# Walk the XML tree, associating every child node with its immediate parent.
parent_map = dict((child, parent) for parent in root.iter() for child in parent)


def ancestor_has_transform(selected_element):
    """Walks up the parent hierarchy for the selected element to determine if there is an associated transformation.

    Arguments:
    selected_element : Some xml element from the trap .svg file.

    Returns:
    Whether the selected element or one of its parents has the transform attribute.
    """
    # As long as the selected element has a parent
    while selected_element in parent_map.keys():
        # See if the selected element has the transform attribute.
        if 'transform' in selected_element.attrib.keys():
            # Stop and report true if we've found a parent with the transform attribute
            return True
        # Otherwise take a step up the tree.
        selected_element = parent_map[selected_element]
    # Finally, if we hit the root element and haven't found the transform attribute, report false.
    return False


def get_matrix(selected_element):
    """Recursive method to return the net transformation matrix of a SVG element

    Arguments:
    selected_element : The element to find the transformation matrix of.

    Returns:
    The net transformation matrix of the selected_element.
    """
    # Start with the identity representing no transformation.
    matrix = np.eye(3)
    # See if the selected matrix has the transform attribute.
    if 'transform' in selected_element.attrib.keys():
        # Get the string representation of the element transformation.
        transform_string = selected_element.attrib['transform']
        # Regular expression representing a float value
        float_regex = r'[-+]?[0-9]*\.?[0-9]+'
        # Create a converted list of floats for all the float values found in the transform_string
        transform_values = [float(val) for val in re.findall(float_regex, transform_string)]

        # The transformation is either a translation, a scaling, or some combination of the two.
        if transform_string.startswith('translate'):
            matrix = np.matrix([[1, 0, transform_values[0]],
                                [0, 1, transform_values[1]],
                                [0, 0, 1]])
        if transform_string.startswith('scale'):
            matrix = np.matrix([[transform_values[0], 0, 0],
                                [0, transform_values[1], 0],
                                [0, 0, 1]])
        if transform_string.startswith('matrix'):
            matrix = np.matrix([[transform_values[0], transform_values[2], transform_values[4]],
                                [transform_values[1], transform_values[3], transform_values[5]],
                                [0, 0, 1]])
    # Check and see if any of the ancestors of our element also have the transformation property.
    if ancestor_has_transform(selected_element):
        # If so, get their combined transformation matrix and multiply it by the one for this element.
        matrix = get_matrix(parent_map[selected_element]) * matrix

    return matrix


class Patch(object):
    """
    Represents a surface patch of an electrode. I think.
    """
    def __init__(self, selected_element):
        self.element = selected_element
        self.id = selected_element.attrib['id']
        self.matched = False

        self.color = '000000'
        if 'style' in selected_element.attrib:
            style = selected_element.attrib['style']
            match = re.search(r'stroke\s*:\s*#(?P<color>\w+)\s*;', style)
            if match:
                self.color = match.group('color')
        self.z = -4.5 if self.color == '0000ff' else 0.0

        if PRINT_EXTRA_OUTPUT:
            print('Creating patch for', self.id)

    def get_centroid(self):
        x, y = 0, 0
        corners = self.get_corners()
        for pos in corners:
            x, y = x + pos[0], y + pos[1]
        return x / len(corners), y / len(corners)

    def distance_to(self, position):
        x, y = self.get_centroid()
        return np.sqrt((x-position[0])**2 + (y-position[1])**2)

    @lru_cache(maxsize=None)
    def get_corners(self):
        transform = get_matrix(parent_map[self.element])
        corners = []
        for p in parse_path(self.element.attrib['d']):
            coordinates = np.matrix([[1, 0, p.start.real], [0, 1, p.start.imag], [0, 0, 1]])
            res = transform*coordinates
            if len(corners):
                dx = corners[-1][0] - res[0, 2]
                dy = corners[-1][1] - res[1, 2]
                if dx*dx + dy*dy < 0.01:
                    continue
            corners.append((res[0, 2], res[1, 2]))
        return corners

    def get_boundaries(self, minx, maxx, miny, maxy):
        corners = self.get_corners()
        prevx, prevy = corners[0]
        inside = minx < prevx < maxx and miny < prevy < maxy

        for pt in corners[1:]:
            x, y = pt
            pt_inside = minx < x < maxx and miny < y < maxy


class Electrode(object):

    all_patches = [Patch(p) for p in root.iterfind('.//g/path')]

    def __init__(self, name):
        self.name = name
        self.patches = [self._get_patch(i)
                        for i in range(len(self.get_label_positions()))]

    @lru_cache(maxsize=None)
    def get_label_positions(self):
        label_positions = []
        for element in root.iterfind('.//tspan'):
            if not element.text == self.name:
                continue
            x = [float(num) for num in element.attrib['x'].split()]  # values are for each letter
            y = [float(num) for num in element.attrib['y'].split()]
            x_avg, y_avg = sum(x)/len(x), sum(y)/len(y)
            final = get_matrix(parent_map[element]) * np.matrix([[1, 0, x_avg], [0, 1, y_avg], [0, 0, 1]])
            label_positions.append((final[0, 2], final[1, 2]))
        return label_positions

    def _get_patch(self, index):
        """Return patches closest to each of the electrodes"""
        delta = float('inf')
        nearest = None
        for patch in Electrode.all_patches:
            if patch.matched:
                continue
            d = patch.distance_to(self.get_label_positions()[index])
            if d < delta:
                delta = d
                nearest = patch
        nearest.matched = True
        return nearest

    all_verts = {}

    def __repr__( self ):
        npoints = [len(p.get_corners()) for p in self.patches]
        s = '{} {}\n'.format(self.name, len(npoints))

        patch_strings = []
        for patch, patch_npoints in zip(self.patches, npoints):
            points = patch.get_corners()
            # for p in points[:-1]:
            # for k, v in self.all_verts.items():
            # if abs( p[0] - v[0] ) < 0.05 and abs( p[1] - v[1] ) < 0.05:
            # print( "OVERLAP: {}-{}".format( k, self.name ) )
            # self.all_verts[ "{}:{}:{}".format( self.name, p[0], p[1] ) ] = p
            points.reverse()
            patch_strings.append('{} {}'.format(patch_npoints, ' '.join( '{} {} {}'.format(
                p[0]*100./9, p[1]*100./9, patch.z) for p in points)))
        return s + '\n'.join(patch_strings)


def main():
    electrode_names = set([elem.text for elem in root.iterfind('.//tspan')])
    electrodes = [Electrode(n) for n in electrode_names]

    xpos, ypos = [], []
    for electrode in electrodes:
        for positions in electrode.get_label_positions():
            xpos.append(positions[0])
            ypos.append(positions[1])

        # print ('Name : ', electrode.name)
        # print ('Labels at : ', electrode.get_label_positions())
        # print ('Patch names : ', [p.id for p in electrode.patches])
        # print ('Patch colors : ', [p.color for p in electrode.patches])
        print(electrode)

    plt.plot(xpos, ypos, 'b.')

    # Plot the labels
    # for x,y,l in zip(xpos, ypos, labels):
    #     plt.text(x,y,l)

    # Plot the shapes
    for patch in Electrode.all_patches:
        # print( patch.get_corners() )
        cent = patch.get_centroid()
        plt.plot(cent[0], cent[1], 'rx')
        # plt.text(cent[0], cent[1], patch.id)
        corners = patch.get_corners()
        x, y = [l[0] for l in corners], [l[1] for l in corners]
        plt.plot(x, y, 'g-')

    plt.show()

if __name__ == "__main__":
    main()
