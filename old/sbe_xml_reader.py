import xml.etree.ElementTree as ET


def sbe_xml_reader_a(file):
    """Function to read .XMLCON file from Seabird.

    Input:
    file: .XMLCON file produced from Seasave 7.

    Output:
    Dictionary of sensors.

    """
    tree = ET.parse(file)
    root = tree.getroot()

    """Pokedex is a dict of {Sensor index numbers from the config:sensor info}
    Assume that sensor index number is the order the sensors have been entered into the file.
    Therefore, it will be Frequency instruments first, then Voltage instruments.
    Call their index number (starting at 0) in order to pull out the info.

    """
    pokedex = {}
    for x in root.iter('Sensor'):
        """Start creating single sensor dictionary."""
        bulbasaur = {}
        for children in x:
            for y in children.iter():
                bulbasaur[y.tag] = float_convert(y.text)
            """Add sensor to big dictionary."""
            pokedex[x.attrib['index']] = bulbasaur
    return pokedex

###### Needs to be looked at ######

def number(string):
    """Helper function to convert strings to int or floats.
    Input: string
    Output: Int, Float or String if non convertable.
    """
    try:
        return int(string)
    except TypeError:
        return float(string)
    except ValueError:
        return string

"""Should be a way to do it with BAFP, need to try later."""
def number_convert(string):
    """Helper function to convert strings to int or floats.
    Input: string
    Output: Int, Float or String if non convertable.
    """
    try:
        return int(string)
    except ValueError:
        try:
            return float(string)
        except TypeError:
            return None
    except ValueError:
        return string

def number_convert_2(string):
    """Use regex to search determine if a string is an int, float, or other.
    Return as int, float or string respectively. Requires "import re".

    Input:
    string: The string to be converted to int, float or string.

    REGEX NOT FULLY TESTED

    """
    print(string)

    #regex patterns
    float_pattern = re.compile('-?[0-9]*\.?[0-9]+|-?[0-9]*\.?[0-9]+e[+-][0-9]*')
    int_pattern = re.compile('-?[0-9]+')

    try:
        if float_pattern.fullmatch(string) is True:
            return float(string)
        if int_pattern.fullmatch(string) is True:
            return int(string)
        return string
    #catch None
    except TypeError:
        return string

###### Needs to be looked at END ######

def float_convert(string):
    try:
        return float(string)
    except:
        return string
