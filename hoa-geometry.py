
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
LOG_EXTRA_OUTPUT = True
EXTRA_OUTPUT_LOG_PATH = 'output/SVGtoExtra.txt'
LOG_ELECTRODES = True
ELECTRODES_LOG_PATH = 'output/SVGtoElectrodes.txt'

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
    Represents the boundary of an electrode.

    The patch object takes in an element from an svg file and parses the
    information into its identity and geometric properties.
    """
    def __init__(self, selected_element, extra_output_log):
        # The id is just an internal name for the element. Used currently for console output.
        self.id = selected_element.attrib['id']
        # The path is the geometric boundary of an electrode.
        self.path = parse_path(selected_element.attrib['d'])
        # All Patches share a common overall scaling, while each has a translation transformation also associated with
        # its relative position in the trap.  The combination of both is represented by the Patch's transform attribute.
        self.transform = get_matrix(parent_map[selected_element])
        # Boolean flag indicating that the Patch has been matched with an electrode identity.
        self.matched = False

        # Extract the style element and grab the stroke information with a regular expression search.
        style = selected_element.attrib['style']
        stroke_regex = r'stroke\s*:\s*#(?P<color>\w+)\s*;'
        stroke_match = re.search(stroke_regex, style)
        # The patch elements are either black ('000000') or blue ('0000ff') which indicates the height
        # of the electrode in the trap. Grab this information from the stroke info and set the height of the patch.
        self.color = stroke_match.group('color')
        # TODO: Figure out why the z value is set to -4.5. The svg is to scale, but I don't know what the scaling is.
        self.z = -4.5 if self.color == '0000ff' else 0.0

        if PRINT_EXTRA_OUTPUT:
            print('Creating patch for', self.id)
        if LOG_EXTRA_OUTPUT:
            extra_output_log.write('Creating patch for ' + self.id + '\n')

    @lru_cache(maxsize=None)  # Allows caching of the corners data to speed up the processing.
    def get_corners(self):
        """Gets the coordinates of the corners of the Patch"""
        # The path specifies a series of line segments, the starting point of each is a corner of the patch
        # traced out by the path.  So we turn each segment starting point into a three-vector and bundle them
        # as a 3x(number of corners) matrix for easy transformation.
        raw_corners = np.transpose(np.array([(segment.start.real, segment.start.imag, 1) for segment in self.path]))
        # The actual corners coordinates are given by the patch transformation applied to the relative coordinates
        # specified by the path geometry data from the previous step.
        corners = self.transform*raw_corners
        return corners

    def get_centroid(self):
        """Finds and returns the centroid of the Patch element.

        The centroid of a patch is simply the average of the coordinate values of the corners.
        """
        corners = self.get_corners()
        # Sum the columns of the corners matrix and divide by the number of columns
        centroid = np.sum(corners, axis=1)/corners[0, :].size
        return centroid

    def distance_to(self, position):
        """Calculates the distance between the centroid of this patch and any position."""
        centroid = self.get_centroid()
        # The third component of the centroid doesn't matter since we're using a 2d projection of the trap.
        distance_vector = centroid[0:2].T-position
        # Standard euclidean norm
        return np.sqrt(np.sum(np.square(distance_vector)))


class Electrode(object):
    """Represents an individual electrode in the trap object.

    An Electrode is an object with a name (as given in the HOA trap manual) and the ability to find the patches
    that represent it within the svg file.  It also knows how to intelligently represent itself when printing
    to the console or a file for further analysis or mesh generation.
    """
    def __init__(self, name, all_patches):
        """From the name and a list of all patches in the SVG file, constructs this electrode."""
        self.name = name
        self.patches = [self._get_patch(i, all_patches) for i in range(len(self.get_label_positions()))]

    @lru_cache(maxsize=None)  # Allows caching of the corners data to speed up the processing.
    def get_label_positions(self):
        """Finds the coordinates of all labels matching the name of this Electrode."""
        label_positions = []
        # Our labels will be carried around in the tspan svg elements.  Loop over these elements and find the
        # ones whose text matches the electrode name.
        for selected_element in root.iterfind('.//tspan'):
            # If we don't match, move on to the next element.
            if not selected_element.text == self.name:
                continue

            # The label has x and y coordinates for each character.
            x = [float(num) for num in selected_element.attrib['x'].split()]
            y = [float(num) for num in selected_element.attrib['y'].split()]
            # Average the character positions to find the center of the label.
            x_avg, y_avg = sum(x)/len(x), sum(y)/len(y)
            # Apply the transformation associated with the element to get its position withing the overall
            # trap geometry.
            final = get_matrix(parent_map[selected_element]) * np.array([[x_avg], [y_avg], [1]])
            label_positions.append([final[0, 0], final[1, 0]])

        return label_positions

    def _get_patch(self, index, all_patches):
        """Return patches closest to each of the electrode labels"""
        # We begin with no distance criteria and no nearest patches.
        delta = float('inf')
        nearest = None
        # Look through all the patches.
        for patch in all_patches:
            # if the patch is already matched to a label, move on.
            if patch.matched:
                continue
            # Find this distance between the current patch and the label we're interested in.
            d = patch.distance_to(self.get_label_positions()[index])
            # If its the smallest distance we've found, make d our distance criteria and the current patch our
            # best candidate for the patch associated with the label of interest.
            if d < delta:
                delta = d
                nearest = patch
        # Once we've found the patch, set it's matched status to true to reduce searching time and return it.
        nearest.matched = True
        return nearest

    def __repr__(self):
        # A list of the number of corners in each patch that makes up the electrode.
        corners_per_patch = [len(patch.get_corners().T) for patch in self.patches]
        # The name followed by the number of patches.
        s = '{} {}\n'.format(self.name, len(self.patches))

        # Place holder for the representation of the patches.
        patch_strings = []
        # TODO: This is another scaling issue that John wrote in but I'm not sure where it comes from.
        scaling = 100./9
        # Loop over all the patches in the electrode
        for patch, number_of_corners in zip(self.patches, corners_per_patch):
            # Grab the corners for this patch
            corners = patch.get_corners().T
            # Now we attach a line with the number of corners followed by a list of their coordinates rescaled and
            # in y, x, z order since the mesh generation code uses a different coordinate geometry.
            # Example for a patch:
            # num_corners corner1_y corner_1_x corner1_z corner2_y corner2_x ... etc.
            patch_strings.append('{} {}'.format(number_of_corners, ' '.join('{} {} {}'.format(
                p[0, 0]*scaling, p[0, 1]*scaling, patch.z) for p in reversed(corners))))
        # Pass back a string with the electrode name followed by the number of patches on the first line,
        # and then a line for each patch as detailed in the loop over the patches.
        return s + '\n'.join(patch_strings)


def main():

    # Open our logging files. This will erase anything in the directory even if logging is set to false.
    extra_output_log = open(EXTRA_OUTPUT_LOG_PATH, 'w')
    electrode_log = open(ELECTRODES_LOG_PATH, 'w')

    # The tspan elements in the svg are all text elements and the only text elements in the svg file are the
    # names of the electrodes. Sort them since using set produces randomly rearranged output files.
    electrode_names = sorted(set([text_element.text for text_element in root.iterfind('.//tspan')]))
    # Grab all of the coordinate information from the svg and create a list of patch objects from each of the
    # path elements.
    all_patches = [Patch(path_element, extra_output_log) for path_element in root.iterfind('.//g/path')]

    # Create electrode objects by pairing the electrode names with all nearest neighbor patches.
    electrodes = [Electrode(name, all_patches) for name in electrode_names]

    # Get the label positions and names for plotting
    label_x_positions, label_y_positions, labels = [], [], []
    for electrode in electrodes:
        for positions in electrode.get_label_positions():
            label_x_positions.append(positions[0])
            label_y_positions.append(positions[1])
            labels.append(electrode.name)

    # Plot the labels
    # As text labels:
    # for x,y,l in zip(label_x_positions, label_y_positions, labels):
    #    plt.text(x,y,l)
    # As blue dots:
    plt.plot(label_x_positions, label_y_positions, 'b.')

    # Plot the patches
    for patch in all_patches:
        # Get the relevant coordinate information for the patch
        corners = np.array(patch.get_corners())
        cent = patch.get_centroid()

        # Plot the patch centroids
        # As the path patch labels (the path from which they were generated):
        # plt.text(cent[0], cent[1], patch.id)
        # As red exes.
        plt.plot(cent[0], cent[1], 'rx')

        # Plot the outlines of the patches as green lines.
        plt.plot(corners[0, :], corners[1, :], 'g-')

        # Optional logging information
        if PRINT_EXTRA_OUTPUT:
            print('Plotting Patch ' + patch.id + ' with corners at : ' + str(corners))
        if LOG_EXTRA_OUTPUT:
            extra_output_log.write('Plotting Patch ' + patch.id + ' with corners at : ' + str(corners) + '\n')

    plt.show()

    # Print information about the patches and electrodes to the console or to a logging file depending on the
    # options selected at the top of the file.
    for electrode in electrodes:
        # Extra output useful for debugging, but unimportant for the mesh generation
        if PRINT_EXTRA_OUTPUT:
            print('Name : ', electrode.name)
            print('Labels at : ', electrode.get_label_positions())
            print('Patch names : ', [p.id for p in electrode.patches])
            print('Patch colors : ', [p.color for p in electrode.patches])
        if LOG_EXTRA_OUTPUT:
            extra_output_log.write('Name : ' + electrode.name + '\n')
            extra_output_log.write('Labels at : ' + str(electrode.get_label_positions()) + '\n')
            extra_output_log.write('Patch names : ' + str([p.id for p in electrode.patches]) + '\n')
            extra_output_log.write('Patch colors : ' + str([p.color for p in electrode.patches]) + '\n')
        # Information used to generate coordinate meshes.
        if PRINT_ELECTRODES:
            print(electrode)
        if LOG_ELECTRODES:
            electrode_log.write(repr(electrode) + '\n')

    extra_output_log.close()
    electrode_log.close()

if __name__ == "__main__":
    main()
