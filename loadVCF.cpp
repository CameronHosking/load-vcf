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

//convert char *s to an unsigned 32bit integer
//len is the number of numeric characters
//s does not require the trailing '\0'
uint32_t atoui32(const char *str, uint8_t len)
{
    uint32_t value = 0;
    switch (len) { // handle up to 10 digits, assume we're 32-bit
            case 10:    value += str[len-10] * 1000000000;
            case  9:    value += str[len- 9] * 100000000;
            case  8:    value += str[len- 8] * 10000000;
            case  7:    value += str[len- 7] * 1000000;
            case  6:    value += str[len- 6] * 100000;
            case  5:    value += str[len- 5] * 10000;
            case  4:    value += str[len- 4] * 1000;
            case  3:    value += str[len- 3] * 100;
            case  2:    value += str[len- 2] * 10;
            case  1:    value += str[len- 1];
        }
    return value - uint32_t(offsets[len]);
}

//convert null terminated char *s to an unsigned 32bit integer
uint32_t atoui32(const char *s)
{
    uint32_t ret = s[0];
    uint8_t len = 1;
    while(s[len])
    {
        ret = ret*10 + s[len++];
    }
    return ret-uint32_t(offsets[len]);
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

// compresses characters as A 00 C 01 G 10 T 11)
//afterwards dest points to the first untouched byte in buffer
//returns wether the compression used long form or short form
//long form uses 4 bits per character and short form uses 2 bits per character
bool writeCompressedString(char * &dest, char * src, size_t length)
{
	
}

struct VariantDetails{
	VariantDetails(const char * text, std::array<size_t,5> positions)
	:pos(atoui32(&(text[positions[0]+1]),positions[1]-positions[0]-1))
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
	
	//writes given alternate allele to buf and returns bytes written
	size_t writeAltToBuf(char * buf, int altNumber)
	{
		const char* alts = getAlt();
		//go to start of nth alternate (counting from 1)
		while(altNumber>1)
		{
			if(*alts==',')
				altNumber--;
			alts++;
		}
		
		int altLength = 0;	
		//add characters that are not \0 or ,
		do{
			buf[altLength] = alts[altLength];
		}while(alts[++altLength]>='A');
		
		return altLength;
	}
	
	//each ploid length is 0 if it is reference
	//second ploid length is 0 if it is homozygous
	//info byte is made up of:
	// number of ploids |homozygous| phased or unphased | is first ploid referance | is second ploid referance |
	// 1 bit            | 1 bit    | 1 bit              | 1 bit                    | 1 bit                     |
	// position in genome | info bit | length of referance | length of ID | first ploid length | second ploid length | ID       | first  ploid seq | second ploid seq 
	// 4 bytes            | 1 byte   | 1 byte              | 1 byte       | 1 byte             | 1 byte              | variable | variable         | variable  
	//returns the amount of bytes written
	size_t writeVariant(char * buf, const char * alleles)
	{
		constexpr uint8_t diploid =     0b10000000;
		constexpr uint8_t homozygous =  0b01000000;
		constexpr uint8_t phased =      0b00100000;
		constexpr uint8_t firstIsRef =  0b00010000;
		constexpr uint8_t secondIsRef = 0b00001000;
		
		//write position in genome
		//force little endian representation (default for x86)
		buf[0] = pos & 0xFF;
		buf[1] = (pos >> 8) & 0xFF;
		buf[2] = (pos >> 16) & 0xFF;
		buf[3] = (pos >> 24) & 0xFF;
		
		const char * c = getId();
		uint8_t idLength = 0;
		//write ID
		while(c[idLength]!='\0')
		{
			buf[9+idLength] = c[idLength];
			idLength++;
		}
		//write ID length
		buf[6] = idLength;
		
		//write reference length
		buf[5] = (uint8_t)strlen(getRef());
		
		
		c = alleles;
		uint8_t infoChar = 0;
		while(*c!='/'&&*c!='|'&&*c!='\0')
		{
			c++;
		}
		
		uint8_t firstPloidLength = 0;
		
		//get first alt number
		size_t firstAltNumber = 0;
		if(c[-1] != '.') 
			firstAltNumber = atoui32(alleles,c-alleles);	
		//write first allele
		if(firstAltNumber == 0)
			infoChar |= firstIsRef;
		else
			firstPloidLength = writeAltToBuf(buf+9+idLength,firstAltNumber);
		buf[7] = firstPloidLength;
		
		//second allele (if diploid)
		uint8_t secondPloidLength = 0;
		if(*c!='\0')//diploid
		{
			infoChar |= diploid;
			if(*c == '|')//phased
				infoChar |= phased;
			
			//get second alt number
			c++;
			size_t secondAltNumber = 0;
			if(*c!='.')
				secondAltNumber = atoui32(c);
			if(secondAltNumber == 0)
				infoChar |= secondIsRef;
			
			if(firstAltNumber==secondAltNumber)//if it is homozygous
				infoChar |= homozygous;
			else if(secondAltNumber != 0)
				secondPloidLength = writeAltToBuf(buf+9+idLength+firstPloidLength,secondAltNumber);
		}
		buf[8] = secondPloidLength;
		
		//write info byte
		buf[4] = infoChar;
		return 9+idLength+firstPloidLength+secondPloidLength;
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

	uint32_t getPos() const
	{
		return pos;
	}
	
	private:
	
	uint32_t pos;
	
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
				break;
			}
		}
	}
	
	void printSampleBinary(std::string id)
	{
		for(auto s: samples)
		{
			if(s.getId()==id)
			{
				size_t buffSize = 1000000;
				char * buf = new char[1000000];
				size_t pos = 0;
				std::ofstream out(id+".bvcf",std::ofstream::out|std::ofstream::binary);
				for(const SampleVariant &v : s.getVariants())
				{
					pos += variantDetails[v.getVariantNumber()].writeVariant(buf+pos,v.getVariantString());
					if(pos > buffSize - 100000)
					{
						out.write(buf,pos);
						pos = 0;
					}
				}
				out.write(buf,pos);
				break;
			}
		}
	}
	
	bool readVCF(bool gzipped = false, const char * filename = "")
	{
		FileReader inFile(1048576,gzipped);
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
	bool gzipped = true;
	vcf.readVCF(gzipped,"minag100k.vcf.gz");
	//vcf.readVCF(gzipped);//read from stdin
	std::vector<std::string> samples({"AA0114-C","AA0050-C","AC0139-C","AJ0038-C","AC0099-C",
	"AB0097-C","AJ0044-C","AB0282-C","AN0039-C","AK0078-C"});
	for(std::string s: samples)
	{
		vcf.printSample(s);
		vcf.printSampleBinary(s);
	}
	return 0;
}
