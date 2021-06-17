#include <iostream>
#include <string>
#include "PartArray.h"
#include "Part.h"
#include "argumentum/include/argumentum/argparse.h"

using namespace std;
using namespace argumentum;

int main(int argc, char* argv[])
{
    auto parser = argumentum::argument_parser{};
    auto params = parser.params();

    std::string filename;

    parser.config().program("distanceAnalyser")
        .description("prints out the distance matrix between \
        all pairs of spins in the system");
    params.add_parameter(filename,"-f","--filename").nargs(1).required().metavar("FILE.mfsys")
        .help("Path to text file with structure of the system. \
            Format is the mfsys file.");

    auto res = parser.parse_args( argc, argv, 1 );

    if ( !res )
      return 1;

    PartArray sys;
    sys.load(filename);

    printf("id");
    for (auto partA:sys.parts){
        printf(";%ld",partA->Id());
    }
    printf("\n");

    for (auto partA:sys.parts){
        printf("%ld",partA->Id());
        for (auto partB:sys.parts){
            printf(";%f",partA->pos.space(partB->pos));
        }
        printf("\n");
    }

    return 0;
}