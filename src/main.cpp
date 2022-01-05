#include "reads_correction_can.h"

int main(int argc, char** argv)
{
    ReadsCorrectionOptions rco;
	int r = parse_arguments(argc, argv, rco);
	if (r) {
		Usage(argv[0]);
		exit(1);
	}
	if (rco.print_usage_info) {
		Usage(argv[0]);
		exit(0);
	}
	
	return reads_correction_can(rco);
}
