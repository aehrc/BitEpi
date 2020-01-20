import argparse

from . import ReadWriteData


def main():
    # Read in the arguments and check for validity
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "-bfile", required=True,
                        help="Argument in the form of: -i <input_file.csv>")
    args = vars(parser.parse_args())
    input_file = args["i"]
    read_write_data = ReadWriteData()
    valid = read_write_data.validate_input_file()
    if valid:
        print("The input file, {}, has been successfully validated."
              .format(input_file))
        read_write_done = read_write_data.read_data_from_csv()
        if read_write_done:
            print(
                "The input file, {}, has been successfully loaded "
                "and the output file has been created successfully.".format(
                    input_file))
        else:
            print("Error has occurred in Read and/or Write of the file.")
    else:
        print("Error found in input file format.")


if __name__ == "__main__":
    main()