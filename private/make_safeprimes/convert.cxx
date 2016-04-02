
#include <icetray/open.h>
#include <icetray/I3Logging.h>

int main (int argc, char const *argv[])
{
	boost::iostreams::filtering_istream ifs;
	boost::iostreams::filtering_ostream ofs;
	
	if (argc < 2)
		log_fatal("Specify one input and one output file!");
	
	I3::dataio::open(ifs, argv[1]);
	if (!ifs.good())
		log_fatal("Could not open input file %s", argv[1]);
	I3::dataio::open(ofs, argv[2], 9);
	if (!ofs.good())
		log_fatal("Could not open output file %s", argv[2]);
	
	// a tag to identify the file
	ofs << "safeprimes_base32";
	
	unsigned i = 0;
	while (!ifs.eof()) {
		int64_t fora;
		ifs >> fora;
		if (ifs.fail()) {
			log_error("Couldn't parse prime at line %u", i+1);
			return 1;
		}
		char c;
		while (!ifs.eof() && ifs.peek() != '\n')
			ifs.read(&c, 1);
		if (fora < std::numeric_limits<uint32_t>::min() || fora > std::numeric_limits<uint32_t>::max())
			break;
		ofs.write(reinterpret_cast<char*>(&fora), sizeof(fora));
		i++;
		if (i % 1000000 == 0)
			log_notice("Wrote %u million primes", i/1000000);
	}
	
	/* code */
	return 0;
}
