#include "grfyreadwrt.h"


namespace Graphene {

	VectorStreamBuf::VectorStreamBuf(std::vector<char>& v) {
		char* begin = v.data();
		size_t size = v.size();
		setbuf(begin, size);
	}
	std::streambuf*  VectorStreamBuf::setbuf(char* s, std::streamsize n)
	{
		setg(s, s, s + n);
		setp(s, s + n);
		return this;
	}
	VectorStreamBuf VectorStream::getVBuf()
	{
		return vbuf;
	}


	//future endian check here too
	//really feel like this should return the value
	



	std::vector<char> file_to_buffer(FILE* fp, size_t to_read, bool locator_change)
	{
		std::vector<char> v(to_read);
		if (to_read == 0)
		{
			fseek(fp, 0, SEEK_END);
			long size = ftell(fp);
			rewind(fp);
			fread(v.data(), sizeof(char), size, fp);
		}
		else {
			fread(v.data(), sizeof(char), to_read, fp);
		}

		if (!locator_change)
		{
			rewind(fp);
		}
		return v;
	}

	inputType determineFType(FILE** fp, char* filename, int index)
	{
		//why is index even needed?
		if (index == 0 || index == 1) {
			if ((fopen_s(fp, filename, "r")) == NULL) {
				//sprintf(buf, "%s_%d", filename, index);
				//if ((*fp = fopen(buf, "r")) == NULL) {
					//fprintf(stderr, "Error, cannot find any file of the form %s\n", filename);
				exit(1);
			}
		}
		std::vector<char> sm = file_to_buffer(*fp, 1, true);
		VectorStream tmp(sm);
		inputType type = *(inputType*)GHRead<int>(tmp);
		return type;

	}

}