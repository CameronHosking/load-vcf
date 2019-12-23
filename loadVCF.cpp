#include <vector>
#include <array>
#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstring>

#include <zlib.h>

static const uint64_t offsets[21] = {0,'0',
        '0'*11ull,
        '0'*111ull,
        '0'*1111ull,
        '0'*11111ull,
        '0'*111111ull,
        '0'*1111111ull,
        '0'*11111111ull,
        '0'*111111111ull,
        '0'*1111111111ull,
        '0'*11111111111ull,
        '0'*111111111111ull,
        '0'*1111111111111ull,
        '0'*11111111111111ull,
        '0'*111111111111111ull,
        '0'*1111111111111111ull,
        '0'*11111111111111111ull,
        '0'*111111111111111111ull,
        '0'*1111111111111111111ull,
        '0'*11111111111111111111ull};

//convert char *s to an unsigned 64bit integer
//len is the number of numeric characters
//s does not require the trailing '\0'
uint64_t atoui64(const char *str, uint8_t len)
{
    size_t value = 0;
    switch (len) { // handle up to 20 digits, assume we're 64-bit
        case 20:    value += str[len-20] * 10000000000000000000ull;
        case 19:    value += str[len-19] * 1000000000000000000ull;
        case 18:    value += str[len-18] * 100000000000000000ull;
        case 17:    value += str[len-17] * 10000000000000000ull;
        case 16:    value += str[len-16] * 1000000000000000ull;
        case 15:    value += str[len-15] * 100000000000000ull;
        case 14:    value += str[len-14] * 10000000000000ull;
        case 13:    value += str[len-13] * 1000000000000ull;
        case 12:    value += str[len-12] * 100000000000ull;
        case 11:    value += str[len-11] * 10000000000ull;
        case 10:    value += str[len-10] * 1000000000ull;
        case  9:    value += str[len- 9] * 100000000ull;
        case  8:    value += str[len- 8] * 10000000ull;
        case  7:    value += str[len- 7] * 1000000ull;
        case  6:    value += str[len- 6] * 100000ull;
        case  5:    value += str[len- 5] * 10000ull;
        case  4:    value += str[len- 4] * 1000ull;
        case  3:    value += str[len- 3] * 100ull;
        case  2:    value += str[len- 2] * 10ull;
        case  1:    value += str[len- 1];
    }
    return value - offsets[len];
}


constexpr size_t SHORT_STRING_LENGTH = 8;
union ShortString
{
	char buf[SHORT_STRING_LENGTH];
	char * buf_ptr;
};

bool setShortString(const char * text, size_t length, ShortString &s)
{
	if(length > SHORT_STRING_LENGTH - 1)
	{
		s.buf_ptr = new char[length+1];
		memcpy(s.buf_ptr,text,length);
		s.buf_ptr[length] = '\0';
		return true;
	}
	else
	{
		memcpy(&(s.buf[0]),text,length);
		s.buf[length] = '\0';
		return false;
	}
}

const char * getShortString(const ShortString &s,bool onHeap)
{
	if(onHeap)
		return s.buf_ptr;
	return &(s.buf[0]);
}

bool copyShortString(ShortString &dest, const ShortString &src, bool srcOnHeap)
{
	const char * s = getShortString(src,srcOnHeap);
	size_t len = strlen(s);
	return setShortString(s,len,dest);
}

struct VariantDetails{
	VariantDetails(const char * text, std::array<size_t,5> positions)
	:pos(atoui64(&(text[positions[0]+1]),positions[1]-positions[0]-1))
	{
		chromOnHeap = setShortString(text,positions[0],chrom);
		idOnHeap = setShortString(&(text[positions[1]+1]),positions[2]-positions[1] - 1, id);
		refOnHeap = setShortString(&(text[positions[2]+1]),positions[3]-positions[2] - 1, ref);
		altOnHeap = setShortString(&(text[positions[3]+1]),positions[4]-positions[3] - 1, alt);
	}
	
	VariantDetails(VariantDetails &&other)
	{
		std::memcpy((void*)this,(void*)(&other),sizeof(VariantDetails));
		other.chromOnHeap = false;
		other.idOnHeap = false;
		other.refOnHeap = false;
		other.altOnHeap = false;
	}
	
	VariantDetails(const VariantDetails &other)
	:pos(other.pos)
	{
		chromOnHeap=copyShortString(chrom,other.chrom,other.chromOnHeap);
		idOnHeap=copyShortString(id,other.id,other.idOnHeap);
		refOnHeap=copyShortString(ref,other.ref,other.refOnHeap);
		altOnHeap=copyShortString(alt,other.alt,other.altOnHeap);
	}
	
	VariantDetails()
	:pos(0)
	{
		chromOnHeap = setShortString("",0,chrom);
		idOnHeap = setShortString("",0, id);
		refOnHeap = setShortString("",0, ref);
		altOnHeap = setShortString("",0, alt);
	}
	
	VariantDetails& operator=(const VariantDetails &other)
	{
		pos = other.pos;
		chromOnHeap=copyShortString(chrom,other.chrom,other.chromOnHeap);
		idOnHeap=copyShortString(id,other.id,other.idOnHeap);
		refOnHeap=copyShortString(ref,other.ref,other.refOnHeap);
		altOnHeap=copyShortString(alt,other.alt,other.altOnHeap);
		return *this;
	}
	
	~VariantDetails()
	{
		if(chromOnHeap)
			delete[] chrom.buf_ptr;
		if(idOnHeap)
			delete[] id.buf_ptr;
		if(refOnHeap)
			delete[] ref.buf_ptr;
		if(altOnHeap)
			delete[] alt.buf_ptr;
	}
		
	const char * getChrom() const
	{
		return getShortString(chrom,chromOnHeap);
	}
	
	const char * getRef() const
	{
		return getShortString(ref,refOnHeap);
	}
	
	const char * getAlt() const
	{
		return getShortString(alt,altOnHeap);
	}
	
	const char * getId() const
	{
		return getShortString(id,idOnHeap);
	}	

	uint64_t getPos() const
	{
		return pos;
	}
	
	private:
	
	uint64_t pos;
	
	bool chromOnHeap;
	bool idOnHeap;
	bool refOnHeap;
	bool altOnHeap;
	
	ShortString chrom;
	ShortString id;
	ShortString ref;
	ShortString alt;
};

	std::ostream& operator<<(std::ostream& o, const VariantDetails &v)
	{
		return (o << v.getChrom() << '\t' << v.getPos() << '\t'
		<< v.getId() << '\t' << v.getRef() << '\t'
		<< v.getAlt());
	}

class FileReader
{
	public:
	const size_t chunkSize;
	char * chunk;
	FILE * inFile;
	gzFile_s * inFile_gz;
	size_t dataToConsume;
	size_t start;
	size_t end;
	bool good;
	bool gzipped;
	public:
	FileReader(size_t chunkSize, bool gzipped = false)
	:chunkSize(chunkSize), chunk(new char[2*chunkSize]),dataToConsume(0),start(0),end(0),good(false),gzipped(gzipped)
	{}
	
	bool open(const char * filename)
	{
		if(gzipped)
		{
			good = ((inFile_gz = gzopen(filename, "rb"))!=NULL);
			return good &= (gzbuffer(inFile_gz,chunkSize)==0);
		}
		return good = ((inFile = fopen(filename, "rb"))!=NULL);
	}
	
	void readFromStdin()
	{
		inFile = stdin;
		good = true;
	}
	
	size_t readIn(char * buf, size_t bytesToRead)
	{
		if(gzipped)
			return gzread(inFile_gz, buf, bytesToRead);
		else
			return fread(buf,1,bytesToRead,inFile);
	}
	
	char readPastDifferentChars(char c1, char c2)
	{
		start = end;
		do{
			while(end < dataToConsume)
			{
				if(chunk[end]==c1||chunk[end]==c2)
				{
					return chunk[end++];
				}
				end++;
			}
			dataToConsume -= start;
			memcpy(chunk,&chunk[start],dataToConsume);
			end = end-start;
			start = 0;
		}while(dataToConsume < chunkSize && (dataToConsume += readIn(chunk+dataToConsume,chunkSize)));
		good = false;
		return '\0';
	}
	template <int N=1>
	bool readPast(char c)
	{
		int num = N;
		start = end;
		do{
			while(end < dataToConsume)
			{
				if(chunk[end++]==c&&((N==1)||(--num==0)))
					return true;
			}
			dataToConsume -= start;
			memcpy(chunk,&chunk[start],dataToConsume);
			end = end-start;
			start = 0;
		}while(dataToConsume < chunkSize && (dataToConsume += readIn(chunk+dataToConsume,chunkSize)));
		good = false;
		return false;
	}
	
	bool readPastAndCheckForCharInBounds(char delim, char lowerBounds, char upperBounds)
	{
		bool flag = false;
		start = end;
		do{
			while(end < dataToConsume)
			{
				if(chunk[end]>=lowerBounds&&chunk[end]<=upperBounds)
					flag = true;				
				if(chunk[end++]==delim)
					return flag;
			}
			dataToConsume -= start;
			memcpy(chunk,&chunk[start],dataToConsume);
			end = end-start;
			start = 0;
		}while(dataToConsume < chunkSize && (dataToConsume += readIn(chunk+dataToConsume,chunkSize)));
		good = false;
		return false;
	}
	
	template <size_t N=1>
	bool skipPast(char c)
	{
		int num = N;
		start = end;
		do{
			while(end < dataToConsume)
			{
				if(chunk[end++]==c&&((N==1)||(--num==0)))
					return true;
			}
			end = 0;
			start = 0;
		}while((dataToConsume = readIn(chunk,chunkSize)));
		good = false;
		return false;
	}
	
	template <size_t N>
	bool markPositions(char c,std::array<size_t,N> &positions)
	{
		start = end;
		int found = 0;
		do{
			while(end < dataToConsume)
			{
				if(chunk[end]==c)
				{
					positions[found++] = end-start;
					end++;
					if(found==N)
						return true;
				}
				end++;
			}
			dataToConsume -= start;
			memcpy(chunk,&chunk[start],dataToConsume);
			end = end-start;
			start = 0;
		}while(dataToConsume < chunkSize && (dataToConsume += readIn(chunk+dataToConsume,chunkSize)));
		good = false;
		return false;
	}
	
	const char* getStartOfRead()
	{
		return &chunk[start];
	}
	
	size_t getCharactersInRead()
	{
		return end - start - 1;
	}
	
	bool isGood(){return good;}
};

class VCF
{	
	class SampleVariant
	{
		uint32_t variantNumber;
		bool onHeap;
		ShortString variant;
		public:
		SampleVariant(uint32_t variantNumber, const char * variantString, size_t variantStringLength)
		:variantNumber(variantNumber)
		{
			onHeap = setShortString(variantString,variantStringLength,variant);
		}
		uint32_t getVariantNumber()const { return variantNumber;}
		const char * getVariantString()const{return getShortString(variant,onHeap);}
	};
	
	class Sample
	{
		std::string id;
		std::vector<SampleVariant> variants;
		public:
		Sample(std::string id)
		:id(id)
		{}
		
		void addVariant(SampleVariant && variant)
		{
			variants.push_back(variant);
		}
		
		const std::string & getId()const{return id;}
		const std::vector<SampleVariant> & getVariants()const{return variants;}
	};
	
	std::vector<VariantDetails> variantDetails;
	std::vector<Sample> samples;
	
	public:
	void printSample(std::string id)
	{
		for(auto s: samples)
		{
			if(s.getId()==id)
			{
				std::ofstream out(id+".vcf");
				out << "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\t" << id <<'\n';
				for(const SampleVariant &v : s.getVariants())
				{
					out << variantDetails[v.getVariantNumber()] << '\t' << v.getVariantString() << '\n';
				}
			}
		}
	}
	bool readVCF(const char * filename = "")
	{
		FileReader inFile(1048576,true);
		if(filename[0]=='\0')
		{
			inFile.readFromStdin();
		}
		else if(!inFile.open(filename))
		{
			std::cout << "failed to open " << filename <<std::endl;
			return false;
		}
		//read in sample IDs line
		char chromTag[] = "#CHROM";
		while(inFile.skipPast('\t'))
		{
			//look at the last 6 characters of each read ending in a positions
			if(memcmp(inFile.getStartOfRead()+inFile.getCharactersInRead()-6,chromTag,6)==0)
			{
				std::cout << "#CHROM " << "found, reading in sample Ids..."<<std::endl;
				inFile.skipPast<8>('\t');
				while(inFile.readPastDifferentChars('\t','\n')=='\t')
				{
					samples.push_back(Sample(std::string(inFile.getStartOfRead(),inFile.getCharactersInRead())));
				}
				//the newline character was read
				if(inFile.isGood())
				{
					samples.push_back(Sample(std::string(inFile.getStartOfRead(),inFile.getCharactersInRead())));
				}
				else//an error occurred (EOF or buffer size exceeded)
				{
					std::cout << "VCF malformed" << std::endl;
					return false;
				}
				break;
			}
		}
		std::cout << "test2"<<std::endl;
		//read in variants for each sample
		while(inFile.isGood())
		{
			//read in variant
			std::array<size_t,5> tabPositions;
			if(inFile.markPositions<5>('\t',tabPositions))
			{
				variantDetails.push_back(VariantDetails(inFile.getStartOfRead(),tabPositions));

				inFile.skipPast<4>('\t');

				for(int i = 0; i < samples.size()-1;++i)
				{
					if(inFile.readPastAndCheckForCharInBounds('\t','1','9'))
					{
						samples[i].addVariant(SampleVariant(variantDetails.size()-1,inFile.getStartOfRead(),inFile.getCharactersInRead()));
					}
				}
				//the last sample is followed by a newline character
				if(inFile.readPastAndCheckForCharInBounds('\n','1','9'))
				{
					samples[samples.size()-1].addVariant(SampleVariant(variantDetails.size()-1,inFile.getStartOfRead(),inFile.getCharactersInRead()));		
				}
			}
		}
		return true;
	}
};
			
								
							

int main()
{
	VCF vcf;
	vcf.readVCF("minag100k.vcf.gz");
	//vcf.readVCF();//read from stdin
	std::vector<std::string> samples({"AA0114-C","AA0050-C","AC0139-C","AJ0038-C","AC0099-C",
	"AB0097-C","AJ0044-C","AB0282-C","AN0039-C","AK0078-C"});
	for(std::string s: samples)
	{
		vcf.printSample(s);
	}
	return 0;
}
