import argparse
from sccore.__init__ import VERSION


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", action="version", version="sccore {}".format(VERSION))
    _args = parser.parse_args()


if __name__ == "__main__":
    main()
