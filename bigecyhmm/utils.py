# Copyright (C) 2024 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Grenoble, France Microcosme
# Univ. Grenoble Alpes, Inria, Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import argparse
import logging
import os
import sys


MIN_VAL = 0
MAX_VAL = 1

logger = logging.getLogger(__name__)


def range_limited_float_type(arg):
    """Type function for argparse - a float within some predefined bounds

    Args:
        arg: argparse argument

    Returns:
        arg: argparse argument
    """
    try:
        f = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a floating point number")
    if f < MIN_VAL or f > MAX_VAL:
        raise argparse.ArgumentTypeError("Argument must be < " + str(MAX_VAL) + " and > " + str(MIN_VAL))
    return f


def limited_integer_type(arg):
    """Type function for argparse - an integer

    Args:
        arg: argparse argument

    Returns:
        arg: argparse argument
    """
    try:
        f = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer number")
    return f


def is_valid_path(filepath):
    """Return True if filepath is valid
    
    Args:
        filepath (str): path to file
    
    Returns:
        bool: True if path exists, False otherwise
    """
    if filepath and not os.access(filepath, os.W_OK):
        try:
            open(filepath, 'w').close()
            os.unlink(filepath)
            return True
        except OSError:
            return False
    else:  # path is accessible
        return True


def is_valid_file(filepath):
    """Return True if filepath exists

    Args:
        filepath (str): path to file

    Returns:
        bool: True if path exists, False otherwise
    """
    try:
        open(filepath, 'r').close()
        return True
    except OSError:
        return False


def is_valid_dir(dirpath):
    """Return True if directory exists or can be created (then create it)
    
    Args:
        dirpath (str): path of directory

    Returns:
        bool: True if dir exists, False otherwise
    """
    if not os.path.isdir(dirpath):
        try:
            os.makedirs(dirpath)
            return True
        except OSError:
            return False
    else:
        return True

def file_or_folder(variable_folder_file):
    """Check if the variable is file or a folder

    Args:
        variable_folder_file (str): path to a file or a folder

    Returns:
        dict: {name of input file: path to input file}
    """
    file_folder_paths = {}

    check_file = False
    if os.path.isfile(variable_folder_file):
        filename = os.path.splitext(os.path.basename(variable_folder_file))[0]
        file_folder_paths[filename] = variable_folder_file
        check_file = True

    check_folder = False
    # For folder, iterate through all files inside the folder.
    if os.path.isdir(variable_folder_file):
        for file_from_folder in os.listdir(variable_folder_file):
            filename = os.path.splitext(os.path.basename(file_from_folder))[0]
            file_folder_paths[filename] = os.path.join(variable_folder_file, file_from_folder)
            check_folder = True

    if check_file is False and check_folder is False:
        logger.critical('ERROR: Wrong input, {0} does not exit'.format(variable_folder_file))
        sys.exit(1)
    return file_folder_paths

