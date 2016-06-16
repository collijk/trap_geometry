
import xml.etree.cElementTree as ET
import re
from matplotlib import pyplot as plt
import numpy as np
from svg.path import parse_path
from functools import lru_cache

it = ET.iterparse('resources/RS1096.svg')
for _, el in it:
    if '}' in el.tag:
        el.tag = el.tag.split('}', 1)[1]  # strip all namespaces
root = it.root
parent_map = dict((c, p) for p in root.getiterator() for c in p)


def ancestor_has_transform(element):
    current = element
    while current in parent_map.keys():
        if 'transform' in current.attrib.keys():
            return True
        current = parent_map[current]
    return False


def get_matrix(element):
    """Recursive method, returns the net transformation matrix of a SVG element"""
    matrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    if 'transform' in element.attrib.keys():
        matrix_text = element.attrib['transform']
        matrix_values = [float(val) for val in re.findall(r'[-+]?[0-9]*\.?[0-9]+', matrix_text)]
        if matrix_text.startswith('matrix'):
            matrix = np.matrix([[matrix_values[0], matrix_values[2], matrix_values[4]],
                                [matrix_values[1], matrix_values[3], matrix_values[5]],
                                [0, 0, 1]])
        if matrix_text.startswith('translate'):
            matrix = np.matrix([[1, 0, matrix_values[0]],
                                [0, 1, matrix_values[1]],
                                [0, 0, 1]])
        if matrix_text.startswith('scale'):
            matrix = np.matrix([[matrix_values[0], 0, 0],
                                [0, matrix_values[1], 0],
                                [0, 0, 1]])
    if ancestor_has_transform(element):
        matrix = get_matrix(parent_map[element]) * matrix
    return matrix


class Patch(object):

    def __init__(self, element):
        self.element = element
        self.id = element.attrib['id']
        self.matched = False

        self.color = '000000'
        if 'style' in element.attrib:
            style = element.attrib[ 'style' ]
            match = re.search( r'stroke\s*:\s*#(?P<color>\w+)\s*;', style )
            if match:
                self.color = match.group( 'color' )
        self.z = -4.5 if self.color == '0000ff' else 0.0

        # print ('Creating patch for', self.id)

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
        for p in parse_path( self.element.attrib['d']):
            coordinates = np.matrix([[1, 0, p.start.real], [0, 1, p.start.imag], [0, 0, 1]])
            res = transform*coordinates
            if len(corners):
                dx = corners[-1][0] - res[0, 2]
                dy = corners[-1][1] - res[1, 2]
                if dx*dx + dy*dy < 0.01:
                    continue
            corners.append((res[0, 2], res[1, 2]))
        return corners

    def get_boundaries( self, minx, maxx, miny, maxy):
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
            if not element.text == self.name: continue
            x = [float(num) for num in element.attrib['x'].split()] # values are for each letter
            y = [float(num) for num in element.attrib['y'].split()]
            x_avg, y_avg = sum(x)/len(x) , sum(y)/len(y)
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