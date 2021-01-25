import tifffile
import numpy as np
import os
from collections import defaultdict
import sys
import re
import h5py

name_regex = re.compile(r"""[\w-]+Pos_(\d+)_(\d+)""")
image_directory = sys.argv[1] 

    
def convert_image_dict_to_h5(h5_dict):     
    for ch_name in h5_dict.keys():
        for pos, image in h5_dict[ch_name].items():
            position_regex = re.findall('\d+', pos)
            major_position = position_regex[0]
            h_name = 'Position_' + major_position + '.h5'
            if os.path.exists(h_name):
                print "Working with existing file: " + h_name
                h = h5py.File(h_name, 'r+')
                if ch_name not in h:
                    group = h.create_group(ch_name)
                else:
                    group = h[ch_name]
                if pos not in h:
                    dataset = group.create_dataset(pos, image.shape, dtype = image.dtype)
                else:
                    dataset = group[pos]
                dataset[...] = image
                h.close()
            else:
                print "Creating New File: " + h_name
                h = h5py.File(h_name, 'w')
                if ch_name not in h:
                    group = h.create_group(ch_name)
                else:
                    group = h[ch_name]
                if pos not in h:
                    dataset = group.create_dataset(pos, image.shape, dtype = image.dtype)
                else:
                    dataset = group[pos]
                dataset[...] = image
                h.close()
    

def get_all_tif_paths(image_directory):
    paths = defaultdict(set)
    for directory, subdirs, filenames in os.walk(image_directory):
        if not filenames:
            continue
        for filename in filenames:
            if not filename.endswith('ome.tif'):
                continue
            paths[directory].add(os.path.join(directory, filename))
    return paths

def get_tif_positions(tif_stack):
    tif_axes = {}
    best_first = 0
    best_second = 0
    for file_path in tif_stack:
        filename = os.path.split(file_path)[1]
        axis_positions = name_regex.search(filename)
        first = int(axis_positions.group(1))
        second = int(axis_positions.group(2))
        best_first = max(first, best_first)
        best_second = max(second, best_second)
        tif_axes[file_path] = (first, second)
    #if best_second > best_first:
        # the second thing is the major axis, so we need to invert them
    #    tif.axes = {file_path: (second, first) for file_path, (first, second) in tif_axes.items()}
    return tif_axes

def create_h5_dictionary(tif_axes):
    h5_dict = defaultdict(lambda: defaultdict(int))
    for file_name, position in tif_axes.items():
        major_position = str(position[0])
        minor_position = str(position[1])
        positions = '(Major, minor) = ({}, {})'.format(major_position, minor_position)
        with tifffile.TiffFile(file_name) as tif:
            image = tif.asarray()
            image_sum = np.sum(image, dtype=np.int, axis=0)
            if 'Cy3' in file_name:
                h5_dict['Cy3'][positions] = image_sum
            if 'Cy5' in file_name:
                h5_dict['Cy5'][positions] = image_sum
    return h5_dict

#def convert_image_dict_to_h5(h5_dict):
#    """Convert nested dictionary to h5 file
#   
#    Parameters
#    ----------
#    name: str
#        Name of the movie; used to name h5 file
#    
#    datadict: dict
#        Dictionary mapping channel names and positions to pixel intensity arrays
#
#    """
#    
#    for ch_name in h5_dict.keys():
#        for pos, image in h5_dict[ch_name].items():
#            position_regex = re.findall('\d+', pos)
#            major_position = position_regex[0]
#            h_name = 'Position_' + major_position + '.h5'
#            if os.path.exists(h_name):
#                print "Working with existing file: " + h_name
#                h = h5py.File('Position_%s.h5' % major_position, 'r+')
#                if ch_name not in h:
#                    group = h.create_group(ch_name)
#                else:
#                    group = h[ch_name]
#                if pos not in h:
#                    dataset = group.create_dataset(pos, image.shape, dtype = image.dtype) 
#                else:
#                    datset = group[pos]
#                dataset[...] = image
#                h.close()
#            else:
#                print "Creating New File: " + h_name
#                h = h5py.File('Position_{}.h5'.format(major_position), 'w')
#                if ch_name not in h:
#                    group = h.create_group(ch_name)
#                else:
#                    group = h[ch_name]
#                if pos not in h:
#                    dataset = group.create_dataset(pos, image.shape, dtype = image.dtype) 
#                else:
#                    datset = group[pos]
#                dataset[...] = image
#def main(image_directory):
if __name__ == '__main__':
    print image_directory
    tif_stack = []
    paths = get_all_tif_paths(image_directory)
    print paths

    for directory, tifs in paths.items():
        tif_stack += list(tifs)
        
    tif_axes = get_tif_positions(tif_stack)
    print tif_axes
    
    h5_dict = create_h5_dictionary(tif_axes)
    
    convert_image_dict_to_h5(h5_dict)
#                h.close()
