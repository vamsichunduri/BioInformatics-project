//David Bierc //Z1737661 //11/25/2018
#include <algorithm> 
#include <tuple>
#include <fstream>//
#include <iostream>
#include <sstream>//
#include <string>//
#include <cstdlib>
#include <locale>
#include <map>
#include <vector>
#include <list>
#include <climits>
using namespace std;
//global variables
const char* fileName = "chr22.hg19.canFam3.sing.maf";
struct STATS 
{
	unsigned int pairMatches = 0;
	unsigned int pairMissmatches = 0;
	unsigned int transitions = 0;
	unsigned int transversions = 0;
	unsigned int gapsTotal = 0;
	std::map <string,unsigned int> matches;
	std::map <unsigned int,unsigned int> gaps;
};

typedef struct STATS Stats;
enum RegionType {UTR3,UTR5,CODINGS,EXONS,INTRONS,PROMOTORS,INTERGENIC};
class Region
{
	int txStart;
	int txEnd;
	vector<tuple<int,int>> exons;
	vector<tuple<int,int>> introns;
	vector<tuple<int,int>> codings;
	vector<tuple<int,int>> utr5;
	vector<tuple<int,int>> utr3;
	vector<tuple<int,int>> promotors;
	vector<tuple<int,int>> intergenic;
	public:
		Region(int start, int end)
		{
			txStart = start;
			txEnd = end;
			intergenic.push_back(tuple<int,int>(start, end));
		}
		//constructor
		Region(int start, int end, int codingStart, int codingEnd, vector<int> starts, vector<int> ends)
		{
			txStart = start;
			txEnd = end;
			exons.reserve(starts.size());
			codings.reserve(starts.size());
			introns.reserve(starts.size()-1);
			for(unsigned int i = 0; i < starts.size(); i++)
			{
				exons.push_back(tuple<int,int>(starts[i],ends[i]));
				//first loop
				if(i == 0)
				{
					codings.push_back(tuple<int,int>(codingStart, ends[i]));
				}
				else if(i == (starts.size()-1))
				{
					codings.push_back(tuple<int,int>(starts[i], codingEnd));
				}
				else
				{
					codings.push_back(tuple<int,int>(starts[i], ends[i]));
				}
				if(i < (starts.size() - 1))
				{
					introns.push_back(tuple<int,int>(ends[i],starts[i+1]));
				}
			}
			promotors.push_back(tuple<int,int>((start - 100), start));
			promotors.push_back(tuple<int,int>((start - 200), start));
			promotors.push_back(tuple<int,int>((start - 400), start));
			promotors.push_back(tuple<int,int>((start - 800), start));
			utr5.push_back(tuple<int,int>(start, codingStart));
			utr3.push_back(tuple<int,int>(codingEnd, end));
		}
		//get methods
		int getStart() const
                {
                        return txStart;
                }
                int getEnd() const
                {
                        return txEnd;
                }
		//methods
		const vector<tuple<int,int>>& getList(RegionType type)
		{
			switch(type)
			{
				case EXONS:
					return exons;
				case CODINGS:
					return codings;
				case INTRONS:
					return introns;
				case UTR3:
					return utr3;
				case UTR5:
					return utr5;
				case PROMOTORS:
					return promotors;
				case INTERGENIC:
					return intergenic;
				default:
					throw; 
			}
		}
		void printRegion()
		{
			cout << "********************************start*****************************" << endl;
			cout << "txStart = " << txStart << " txEnd = " << txEnd << endl;
			cout << "Exons:";
			for(unsigned int i = 0; i < exons.size(); i++)
			{
				cout << " " << get<0>(exons[i]) << "," << get<1>(exons[i]);
			}
			cout << endl;
			cout << "Codings:";
                        for(unsigned int i = 0; i < codings.size(); i++)
                        {
                                cout << " " << get<0>(codings[i]) << "," << get<1>(codings[i]);
                        }
                        cout << endl;
			cout << "Introns:";
                        for(unsigned int i = 0; i < introns.size(); i++)
                        {
                                cout << " " << get<0>(introns[i]) << "," << get<1>(introns[i]);
                        }
                        cout << endl;
			cout << "Utr3:";
                        for(unsigned int i = 0; i < utr3.size(); i++)
                        {
                                cout << " " << get<0>(utr3[i]) << "," << get<1>(utr3[i]);
                        }
                        cout << endl;
			cout << "Utr5:";
                        for(unsigned int i = 0; i < utr5.size(); i++)
                        {
                                cout << " " << get<0>(utr5[i]) << "," << get<1>(utr5[i]);
                        }
                        cout << endl;
			cout << "Promoters:";
                        for(unsigned int i = 0; i < promotors.size(); i++)
                        {
                                cout << " " << get<0>(promotors[i]) << "," << get<1>(promotors[i]);
                        }
                        cout << endl;
			cout << "Intergenic:";
                        for(unsigned int i = 0; i < intergenic.size(); i++)
                        {
                                cout << " " << get<0>(intergenic[i]) << "," << get<1>(intergenic[i]);
                        }
                        cout << endl;
			cout << "********************************End*******************************" << endl;
		}
};

vector<Region*> Regions;
void createIntergincRegions(int intergenicStart, int intergenicEnd)
{
	struct 
	{
		bool operator()(const Region* a, const Region* b)
		{
			return a->getStart() < b->getStart();
		}
	}RegionsCompare;

        //sort regions by region
	std::sort(Regions.begin(), Regions.end(),RegionsCompare); 

	vector<Region*> regionList2;
        //add new regions
	if (intergenicStart < Regions[0]->getStart())
            	regionList2.push_back(new Region(intergenicStart,Regions[0]->getStart()));
	for(unsigned int i = 0; i < (Regions.size() - 1); i++)
	{
		if ((Regions[i]->getEnd()) < Regions[i+1]->getStart())
			regionList2.push_back(new Region(Regions[i]->getEnd(), Regions[i+1]->getStart()));
	}
	if(intergenicEnd > Regions[Regions.size()-1]->getEnd())
               regionList2.push_back(new Region(intergenicEnd,Regions[Regions.size()-1]->getEnd()));
	//cobine the new array
	for(unsigned int i = 0; i < regionList2.size(); i++)
	{
		Regions.push_back(regionList2[i]);
	}
}

vector<tuple<int,int>> getRegions(int start, int end , RegionType type)
{
	vector<tuple<int,int>>  overlapping;
	for(unsigned int i = 0; i < Regions.size(); i++)
	{
		//check
		if( Regions[i]->getStart() > end)
			continue;
		if( Regions[i]->getEnd() < start)
			continue;
		const vector<tuple<int,int>>& list = Regions[i]->getList(type);
		for(unsigned int j = 0; j < list.size();j++)
		{
			if(get<1>(list[j]) < start)
                        	continue;
			if(get<0>(list[j]) > end)
                        	continue;
			overlapping.push_back(list[j]);
		}
	}
	return overlapping;
}
//function PrintGaps used to print out map of gaps in readable manner
void printGaps(Stats &s)
{
	//iterate threw the map of gaps one by one
	for(map<unsigned int, unsigned int>::const_iterator it = s.gaps.begin();it != s.gaps.end(); ++it)
        {
                cout << "gaps of length: " << it->first << " Count: " << it->second << " Gap Frequency(Gap length/gaps Total): " << double(it->second)/(s.gapsTotal) << endl;
        }
	cout << "Total number of gaps = " << s.gapsTotal << endl; 
}
//function printMatches used to print out map of matches in readable manner
void printMatches(Stats &s)
{
	//iterate threw the map of matches one by one
	for(map<string, unsigned int>::const_iterator it = s.matches.begin();it != s.matches.end(); ++it)
	{
   		cout << "pair: " << it->first << " Count: " << it->second << endl;
		//get all matches
		if((it->first == "AA") ||(it->first == "CC") ||(it->first == "GG") ||(it->first == "TT"))
		{
			s.pairMatches += it->second;
		}
		//get all missmatches
		if((it->first == "AC") || (it->first == "AG") || (it->first == "AT") || (it->first == "CA") || (it->first == "GA") || (it->first == "TA") || (it->first == "GC") || (it->first == "TC") || (it->first == "TG") || (it->first == "CG") || (it->first == "CT") || (it->first == "GT"))
		{
			s.pairMissmatches += it->second;
		}
		//get transitions
		if((it->first == "AG") || (it->first == "GA") ||  (it->first == "CT") || (it->first == "TC"))
		{
			s.transitions += it->second;
		}
		//get Transversions
		if( (it->first == "AC") || (it->first == "CA") || (it->first == "AT") || (it->first == "TA") || (it->first == "GC") || (it->first == "CG") || (it->first == "GT") || (it->first == "TG"))
		{
			s.transversions += it->second;
		}
	}
	//calculate 
	double subsitutionRate = double(s.pairMissmatches) / (s.pairMissmatches + s.pairMatches);
	double titv = double(s.transitions) / s.transversions;
	double gapRate = double(s.gapsTotal) / (s.pairMissmatches + s.pairMatches + s.gapsTotal);
	//display all the counts and calculations
	cout << "Matches: " << s.pairMatches << endl;
	cout << "Mismatches: " << s.pairMissmatches << endl;
	cout << "Transitions: " << s.transitions << endl;
	cout << "Transversions: " << s.transversions << endl;
	cout << "Substitution Rate: " << subsitutionRate << endl;
	cout << "TiTv: " << titv << endl;
	cout << "Gap Rate: " << gapRate << endl;
}
std::vector<string> stringSplit(string &line)
{
        std::vector<string> result;
	stringstream check1(line);
	string intermediate;
	while(getline(check1, intermediate, ' ')) 
    	{
		//don't include the empty strings that are created by the weird multiple spaces
		if(intermediate != "")
		{
        		result.push_back(intermediate);
		}
    	}
        return result;
}
std::vector<string> stringSplit2(string &line)
{
        std::vector<string> result;
        stringstream check1(line);
        string intermediate;
        while(getline(check1, intermediate, ','))
        {
                //don't include the empty strings that are created by the weird multiple spaces
                if(intermediate != "")
                {
                        result.push_back(intermediate);
                }
        }
        return result;
}
std::vector<string> RemoveGapsInSpiecies1AndIndexofSpecies2(string &line1, string &line2)
{
	int index = 0;
	std::vector<string> result;
	//erase the gaps from the human genome
	while((index = line1.find("-")) != std::string::npos)
	{
		line1.erase(index,1);
		line2.erase(index,1);
	}
	result.push_back(line1);
	result.push_back(line2); 
	return result;
}
//fuunction is used specifically to put the comparision of pairs into map data containers
void compare(string &line1, string &line2, Stats &s)
{
	//local varibales
	unsigned int i = 0;
	unsigned int gapLength1 = 0;
	unsigned int gapLength2 = 0;
	bool isPerviousGap1 = false;
	bool isPerviousGap2 = false;
	//since both lines are cheacked to be equal I just picked one length and iterated char by char
	while(i < line1.length())
	{
		//grab the 2 characters
		char s1 = line1[i];
		char s2 = line2[i];
		//make sure they are uppercase
		s1 = toupper(s1);
		s2 = toupper(s2);
		//build the key
		string key;
                key = s1;
                key += s2;
		//check for end of gap in line 1
		if( isPerviousGap1 && (s1 != '-'))
		{
			//gap has come to an end
			unsigned int key = gapLength1;
			//check for key
			if(s.gaps.count(key) > 0)
                	{
                        	s.gaps.at(key) = s.gaps.at(key) + 1;
				s.gapsTotal++;
                	}
                	else
			{
				//generate new key in data container
				s.gaps.insert(pair<unsigned int,unsigned int>(key, 1));
				s.gapsTotal++;
			}
			//reset
			isPerviousGap1 = false;
			gapLength1 = 0;
		}
		//check of end of gap line 2
		if( isPerviousGap2 && (s2 != '-'))
                {
                        //gap has come to an end
                        unsigned int key = gapLength2;
			//check for key
                        if(s.gaps.count(key) > 0)
                        {
                                s.gaps.at(key) = s.gaps.at(key) + 1;
				s.gapsTotal++;
                        }
                        else
                        {
				//generate new key in data container
                                s.gaps.insert(pair<unsigned int,unsigned int>(key, 1));
				s.gapsTotal++;
                        }
			//reset
                        isPerviousGap2 = false;
			gapLength2 = 0;
                }

		//check if gap case in line 1
		if(s1 == '-')
		{
			gapLength1++;
			isPerviousGap1 = true;
		}
		//check for gap space in line 2
		if(s2 == '-')
		{
			gapLength2++;
			isPerviousGap2 = true;
		}
		
		//generates matches
		if(s.matches.count(key) > 0)
		{
			s.matches.at(key) = s.matches.at(key) + 1;
		}
		else 
		{
			//generate new key in data container
			s.matches.insert(pair<string,unsigned int>(key, 1));
		}
		i++;
	}
}
void readtxtFile()
{
	string line;
	ifstream inFile("hg19.knownGene.txt");
        if(!inFile)
        {
                cerr << "Unable to open file";
        }
        //first info ignored
        getline(inFile,line);
	//variables for storage
	string name;
	string chrom;
	string strand;
	int txStart;
	int txEnd;
	int cdsStart;
	int cdsEnd;
	int exonCount;
	string exonStart;
	string exonEnd;
	unsigned int startIntergenic = 0;
	unsigned int endIntergenic = INT_MAX;
	int start22 = INT_MAX;
	int end22 = 0;
	while(getline(inFile,line))
	{
		istringstream iss(line);
		if(!(iss >> name >> chrom >> strand >> txStart >> txEnd >> cdsStart >> cdsEnd >> exonCount >> exonStart >> exonEnd))
		{
			cout << "Error reading Txt file(first line)." << endl; 
			break;//Error
		}
		//look for the right chromosome
		if(chrom == "chr22")
		{
			if(txStart < start22)
			{
				start22 = txStart;
			}
			if(txEnd > end22)
			{
				end22 = txEnd;
			}
			//split the commas from the string to get all the integers
			vector<string> exonStartTokens = stringSplit2(exonStart);
                	vector<string> exonEndTokens = stringSplit2(exonEnd);
			
			//convert string vector to int vector
			vector<int> exonStart;
                        vector<int> exonEnd;
			for(unsigned int i = 0; i < exonStartTokens.size(); i++)
			{
				//convert and push into new array
				exonStart.push_back(stoi(exonStartTokens[i]));
				exonEnd.push_back(stoi(exonEndTokens[i]));
			}
			//put it into region vector
			Regions.push_back(new Region(txStart, txEnd, cdsStart, cdsEnd, exonStart, exonEnd));
		}
	}
	inFile.close();
	ifstream inFile2("hg19.knownGene.txt");
        if(!inFile2)
        {
                cerr << "Unable to open file";
        }
        //first info ignored
        getline(inFile2,line);
	//restart search
	//inFile.seekg(0, inFile.beg);
	while(getline(inFile2,line))
        {
                istringstream iss(line);
                if(!(iss >> name >> chrom >> strand >> txStart >> txEnd >> cdsStart >> cdsEnd >> exonCount >> exonStart >> exonEnd))
		{
			cout << "Error reading Txt file(second line)." << endl;
                        break;//Error
		}
                //look for the right chromosome
                if(chrom != "chr22")
                {
			if(txEnd > startIntergenic && txEnd <= start22)
				startIntergenic = txEnd;
			if(txStart < endIntergenic && txStart >= end22)
				endIntergenic = txStart;
		}

	}

	inFile2.close();
	vector<Region*> Intergenic;
	for(unsigned int i =0; i < Regions.size(); i++)
	{
		unsigned int txStart = Regions[i]->getStart();
		unsigned int txEnd = Regions[i]->getEnd();
		if( i == 0) //left most region
		{
			if(startIntergenic < txStart)//check left side for intergenic gap
			{
				Intergenic.push_back(new Region(startIntergenic, (txStart - 1)));// create left most region if valid
			}
		}
		
		if(i == (Regions.size()-1)) //right most region
		{
			if(endIntergenic > txEnd)//check right side for intergenic gap
			{
				Intergenic.push_back(new Region((txEnd + 1), endIntergenic)); // create right most region if vaild 
			}
		}
		else //normal case in between already discovered regions
		{
			unsigned int txStartNext = Regions[i+1]->getStart(); //grab next region
			unsigned int txEndNext = Regions[i+1]->getEnd();
			Intergenic.push_back(new Region((txEnd + 1),(txStartNext-1)));
		}
	}
	for(unsigned int i =0; i < Intergenic.size(); i++)
	{
		Regions.push_back(Intergenic[i]);
	}
	//sort it for debugging  
	struct
        {
                bool operator()(const Region* a, const Region* b)
                {
                        return a->getStart() < b->getStart();
                }
        }RegionsCompare;

        //sort regions by region
        std::sort(Regions.begin(), Regions.end(),RegionsCompare);

	cout << "Done reading txt file" << endl;
}
//arguments argv used to get file names
int main(int argc, char **argv)
{
	readtxtFile();
	//print regions for debugging table
	//for(unsigned int i =0 ; i < 7; i++)
	//{
	//	Regions[i]->printRegion();
	//}

	string line;

	ifstream inFile(fileName);
	if(!inFile)
	{
		cerr << "Unable to open file";
	}
	//first info ignored
	getline(inFile,line);

	int pairCount = 1;
	char ignor = ' ';
	string speciesChrom1 = "";
	string speciesChrom2 = "";
	unsigned int sequenceOneStartPosition = 0;
	unsigned int sequenceTwoStartPosition = 0;
	unsigned int sequenceOneLength = 0;
        unsigned int sequenceTwoLength = 0;
	char strand = ' '; //also ignored
	unsigned int chromosomeSize1 = 0;
	unsigned int chromosomeSize2 = 0;
	string alignmentSequence1 = "";
	string alignmentSequence2 = "";

	Stats ExonsStats;
	Stats Utr3Stats;
	Stats Utr5Stats;
	Stats CodingsStats;
	Stats IntronsStats;
	Stats IntergenicStats;
	Stats Promoters100Stats;
	Stats Promoters200Stats;
	Stats Promoters400Stats;
	Stats Promoters800Stats;

	while(getline(inFile,line))
	{
		//reset variables
		ignor = ' ';
        	speciesChrom1 = "";
        	speciesChrom2 = "";
        	sequenceOneStartPosition = 0;
        	sequenceTwoStartPosition = 0;
        	sequenceOneLength = 0;
        	sequenceTwoLength = 0;
        	strand = ' '; //also ignored
        	chromosomeSize1 = 0;
        	chromosomeSize2 = 0;
        	alignmentSequence1 = "";
        	alignmentSequence2 = "";

		if(line[0] == 's')
		{
			//cout << line << endl;
			istringstream iss(line);
                	if(!(iss >> ignor >> speciesChrom1 >> sequenceOneStartPosition >> sequenceOneLength >> strand >> chromosomeSize1 >> alignmentSequence1))
			{
				cout << "Error reading maf file(first line)." << endl;
				break;//Error
			}

			//grab and test second sequence
			getline(inFile,line);
			istringstream iss2(line);
                        if(!(iss2 >> ignor >> speciesChrom2 >> sequenceTwoStartPosition >> sequenceTwoLength >> strand >> chromosomeSize2 >> alignmentSequence2))
			{
				cout << "Error reading maf file(second line)." << endl;
				break;//Error
			}
			//to upper
			transform(alignmentSequence1.begin(), alignmentSequence1.end(), alignmentSequence1.begin(), ::toupper);
			transform(alignmentSequence2.begin(), alignmentSequence2.end(), alignmentSequence2.begin(), ::toupper);

			//align sequence to human length so that positions can be match knownhumangene.txt 
			vector<string> alignmentsWithGapsRemoved = RemoveGapsInSpiecies1AndIndexofSpecies2(alignmentSequence1, alignmentSequence2);

			unsigned int sequenceOneEndPosition = sequenceOneStartPosition + sequenceOneLength;
			//look up for Exons with in this pair of .maf file**********************************************************************************
			vector<tuple<int,int>> exonRegions = getRegions(sequenceOneStartPosition, sequenceOneEndPosition, EXONS);
			if(exonRegions.size() == 0)
			{
				//cout << "found nothing overlapping for ....Exons.... in this pair in maf file." << endl;
			}
			else
			{
				//do calculations for each exon
				for(unsigned int i = 0; i < exonRegions.size(); i++)
				{
					//get the start and end
					unsigned int regionStart = get<0>(exonRegions[i]);
					unsigned int regionEnd = get<1>(exonRegions[i]);

					unsigned int start = 0;
					unsigned int end = 0;
					//left side overlapping exon region
					if(regionStart > sequenceOneStartPosition)
					{
						start = regionStart - sequenceOneStartPosition;
					}
					else
					{
						start = 0;
					}
					//right side overlapping exon region
					if(regionEnd > sequenceOneEndPosition)
					{
						end = sequenceOneLength;
					}
					else
					{
						end = regionEnd - sequenceOneStartPosition;
					}
					//cut the region out
					string alignmentSeq1Sub = alignmentsWithGapsRemoved[0].substr(start, end);
					string alignmentSeq2Sub = alignmentsWithGapsRemoved[1].substr(start, end);
					//run calculations on sub strings
					compare(alignmentSeq1Sub, alignmentSeq2Sub, ExonsStats);
				}
			}
			//End of look up for Exons *********************************************************************************************************************
			//look up for introns with in this pair of .maf file**********************************************************************************
                        vector<tuple<int,int>> intronRegions = getRegions(sequenceOneStartPosition, sequenceOneEndPosition, INTRONS);
                        if(intronRegions.size() == 0)
                        {
                                //cout << "found nothing overlapping for ....Introns.... in this pair in maf file." << endl;
                        }
                        else
                        {
                                //do calculations for each exon
                                for(unsigned int i = 0; i < intronRegions.size(); i++)
                                {
                                        //get the start and end
                                        unsigned int regionStart = get<0>(intronRegions[i]);
                                        unsigned int regionEnd = get<1>(intronRegions[i]);

                                        unsigned int start = 0;
                                        unsigned int end = 0;
                                        //left side overlapping exon region
                                        if(regionStart > sequenceOneStartPosition)
                                        {
                                                start = regionStart - sequenceOneStartPosition;
                                        }
                                        else
                                        {
                                                start = 0;
                                        }
                                        //right side overlapping exon region
                                        if(regionEnd > sequenceOneEndPosition)
                                        {
                                                end = sequenceOneLength;
                                        }
                                        else
                                        {
                                                end = regionEnd - sequenceOneStartPosition;
                                        }
                                        //cut the region out
                                        string alignmentSeq1Sub = alignmentsWithGapsRemoved[0].substr(start, end);
                                        string alignmentSeq2Sub = alignmentsWithGapsRemoved[1].substr(start, end);
                                        //run calculations on sub strings
                                        compare(alignmentSeq1Sub, alignmentSeq2Sub, IntronsStats);
                                }
                        }
                        //End of look up for Introns *********************************************************************************************************************
			//look up for Codings with in this pair of .maf file**********************************************************************************
                        vector<tuple<int,int>> codingRegions = getRegions(sequenceOneStartPosition, sequenceOneEndPosition, CODINGS);
                        if(codingRegions.size() == 0)
                        {
                                //cout << "found nothing overlapping for ....Codings.... in this pair in maf file." << endl;
                        }
                        else
                        {
                                //do calculations for each exon
                                for(unsigned int i = 0; i < codingRegions.size(); i++)
                                {
                                        //get the start and end
                                        unsigned int regionStart = get<0>(codingRegions[i]);
                                        unsigned int regionEnd = get<1>(codingRegions[i]);

                                        unsigned int start = 0;
                                        unsigned int end = 0;
                                        //left side overlapping exon region
                                        if(regionStart > sequenceOneStartPosition)
                                        {
                                                start = regionStart - sequenceOneStartPosition;
                                        }
                                        else
                                        {
                                                start = 0;
                                        }
                                        //right side overlapping exon region
                                        if(regionEnd > sequenceOneEndPosition)
                                        {
                                                end = sequenceOneLength;
                                        }
                                        else
                                        {
                                                end = regionEnd - sequenceOneStartPosition;
                                        }
                                        //cut the region out
                                        string alignmentSeq1Sub = alignmentsWithGapsRemoved[0].substr(start, end);
                                        string alignmentSeq2Sub = alignmentsWithGapsRemoved[1].substr(start, end);
                                        //run calculations on sub strings
                                        compare(alignmentSeq1Sub, alignmentSeq2Sub, CodingsStats);
                                }
                        }
                        //End of look up for codings *********************************************************************************************************************
			//look up for utr3 with in this pair of .maf file**********************************************************************************
                        vector<tuple<int,int>> utr3Regions = getRegions(sequenceOneStartPosition, sequenceOneEndPosition, UTR3);
                        if(utr3Regions.size() == 0)
                        {
                                //cout << "found nothing overlapping for ....UTR3.... in this pair in maf file." << endl;
                        }
                        else
                        {
                                //do calculations for each exon
                                for(unsigned int i = 0; i < utr3Regions.size(); i++)
                                {
                                        //get the start and end
                                        unsigned int regionStart = get<0>(utr3Regions[i]);
                                        unsigned int regionEnd = get<1>(utr3Regions[i]);

                                        unsigned int start = 0;
                                        unsigned int end = 0;
                                        //left side overlapping exon region
                                        if(regionStart > sequenceOneStartPosition)
                                        {
                                                start = regionStart - sequenceOneStartPosition;
                                        }
                                        else
                                        {
                                                start = 0;
                                        }
                                        //right side overlapping exon region
                                        if(regionEnd > sequenceOneEndPosition)
                                        {
                                                end = sequenceOneLength;
                                        }
                                        else
                                        {
                                                end = regionEnd - sequenceOneStartPosition;
                                        }
                                        //cut the region out
                                        string alignmentSeq1Sub = alignmentsWithGapsRemoved[0].substr(start, end);
                                        string alignmentSeq2Sub = alignmentsWithGapsRemoved[1].substr(start, end);
                                        //run calculations on sub strings
                                        compare(alignmentSeq1Sub, alignmentSeq2Sub, Utr3Stats);
                                }
                        }
                        //End of look up for utr3 *********************************************************************************************************************
			//look up for Codings with in this pair of .maf file**********************************************************************************
                        vector<tuple<int,int>> utr5Regions = getRegions(sequenceOneStartPosition, sequenceOneEndPosition, UTR5);
                        if(utr5Regions.size() == 0)
                        {
                                //cout << "found nothing overlapping for ....UTR5.... in this pair in maf file." << endl;
                        }
                        else
                        {
                                //do calculations for each exon
                                for(unsigned int i = 0; i < utr5Regions.size(); i++)
                                {
                                        //get the start and end
                                        unsigned int regionStart = get<0>(utr5Regions[i]);
                                        unsigned int regionEnd = get<1>(utr5Regions[i]);

                                        unsigned int start = 0;
                                        unsigned int end = 0;
                                        //left side overlapping exon region
                                        if(regionStart > sequenceOneStartPosition)
                                        {
                                                start = regionStart - sequenceOneStartPosition;
                                        }
                                        else
                                        {
                                                start = 0;
                                        }
                                        //right side overlapping exon region
                                        if(regionEnd > sequenceOneEndPosition)
                                        {
                                                end = sequenceOneLength;
                                        }
                                        else
                                        {
                                                end = regionEnd - sequenceOneStartPosition;
                                        }
                                        //cut the region out
                                        string alignmentSeq1Sub = alignmentsWithGapsRemoved[0].substr(start, end);
                                        string alignmentSeq2Sub = alignmentsWithGapsRemoved[1].substr(start, end);
                                        //run calculations on sub strings
                                        compare(alignmentSeq1Sub, alignmentSeq2Sub, Utr5Stats);
                                }
                        }
                        //End of look up for codings *********************************************************************************************************************
			//look up for promotors with in this pair of .maf file**********************************************************************************
                        vector<tuple<int,int>> promotorRegions = getRegions(sequenceOneStartPosition, sequenceOneEndPosition, PROMOTORS);
                        if(promotorRegions.size() == 0)
                        {
                                //cout << "found nothing overlapping for ....Interginic.... in this pair in maf file." << endl;
                        }
                        else
                        {
				//ptomotor 100 start*********************************************
				//get the start and end
                                unsigned int regionStart = get<0>(promotorRegions[0]);
                                unsigned int regionEnd = get<1>(promotorRegions[0]);

                                unsigned int start = 0;
                                unsigned int end = 0;
                                //left side overlapping exon region
                                if(regionStart > sequenceOneStartPosition)
                                {
                                       start = regionStart - sequenceOneStartPosition;
                                }
                                else
                                {
                                       start = 0;
                                }
                                //right side overlapping exon region
                                if(regionEnd > sequenceOneEndPosition)
                                {
                                    end = sequenceOneLength;
                                }
                                else
                                {
                                       end = regionEnd - sequenceOneStartPosition;
                                }
                                //cut the region out
                                string alignmentSeq1Sub = alignmentsWithGapsRemoved[0].substr(start, end);
                                string alignmentSeq2Sub = alignmentsWithGapsRemoved[1].substr(start, end);
                                //run calculations on sub strings
                                compare(alignmentSeq1Sub, alignmentSeq2Sub, Promoters100Stats);
				//promotor 100 end*********************************************
				
				//ptomotor 200 start*********************************************
                                //get the start and end
                                regionStart = get<0>(promotorRegions[1]);
                                regionEnd = get<1>(promotorRegions[1]);

                                start = 0;
                                end = 0;
                                //left side overlapping exon region
                                if(regionStart > sequenceOneStartPosition)
                                {
                                       start = regionStart - sequenceOneStartPosition;
                                }
                                else
                                {
                                       start = 0;
                                }
                                //right side overlapping exon region
                                if(regionEnd > sequenceOneEndPosition)
                                {
                                    end = sequenceOneLength;
                                }
                                else
                                {
                                       end = regionEnd - sequenceOneStartPosition;
                                }
                                //cut the region out
                                alignmentSeq1Sub = alignmentsWithGapsRemoved[0].substr(start, end);
                                alignmentSeq2Sub = alignmentsWithGapsRemoved[1].substr(start, end);
                                //run calculations on sub strings
                                compare(alignmentSeq1Sub, alignmentSeq2Sub, Promoters200Stats);
                                //promotor 200 end*********************************************
				
				//ptomotor 400 start*********************************************
                                //get the start and end
                                regionStart = get<0>(promotorRegions[2]);
                                regionEnd = get<1>(promotorRegions[2]);

                                start = 0;
                                end = 0;
                                //left side overlapping exon region
                                if(regionStart > sequenceOneStartPosition)
                                {
                                       start = regionStart - sequenceOneStartPosition;
                                }
                                else
                                {
                                       start = 0;
                                }
                                //right side overlapping exon region
                                if(regionEnd > sequenceOneEndPosition)
                                {
                                    end = sequenceOneLength;
                                }
                                else
                                {
                                       end = regionEnd - sequenceOneStartPosition;
                                }
                                //cut the region out
                                alignmentSeq1Sub = alignmentsWithGapsRemoved[0].substr(start, end);
                                alignmentSeq2Sub = alignmentsWithGapsRemoved[1].substr(start, end);
                                //run calculations on sub strings
                                compare(alignmentSeq1Sub, alignmentSeq2Sub, Promoters400Stats);
                                //promotor 400 end*********************************************

				//ptomotor 800 start*********************************************
                                //get the start and end
                                regionStart = get<0>(promotorRegions[3]);
                                regionEnd = get<1>(promotorRegions[3]);

                                start = 0;
                                end = 0;
                                //left side overlapping exon region
                                if(regionStart > sequenceOneStartPosition)
                                {
                                       start = regionStart - sequenceOneStartPosition;
                                }
                                else
                                {
                                       start = 0;
                                }
                                //right side overlapping exon region
                                if(regionEnd > sequenceOneEndPosition)
                                {
                                    end = sequenceOneLength;
                                }
                                else
                                {
                                       end = regionEnd - sequenceOneStartPosition;
                                }
                                //cut the region out
                                alignmentSeq1Sub = alignmentsWithGapsRemoved[0].substr(start, end);
                                alignmentSeq2Sub = alignmentsWithGapsRemoved[1].substr(start, end);
                                //run calculations on sub strings
                                compare(alignmentSeq1Sub, alignmentSeq2Sub, Promoters800Stats);
                                //promotor 800 end*********************************************
                        }
                        //End of look up for intergenics *********************************************************************************************************************

			//look up for intergenics with in this pair of .maf file**********************************************************************************
                        vector<tuple<int,int>> intergenicRegions = getRegions(sequenceOneStartPosition, sequenceOneEndPosition, INTERGENIC);
                        if(intergenicRegions.size() == 0)
                        {
                                //cout << "found nothing overlapping for ....Interginic.... in this pair in maf file." << endl;
                        }
                        else
                        {
                                //do calculations for each exon
                                for(unsigned int i = 0; i < intergenicRegions.size(); i++)
                                {
                                        //get the start and end
                                        unsigned int regionStart = get<0>(intergenicRegions[i]);
                                        unsigned int regionEnd = get<1>(intergenicRegions[i]);

                                        unsigned int start = 0;
                                        unsigned int end = 0;
                                        //left side overlapping exon region
                                        if(regionStart > sequenceOneStartPosition)
                                        {
                                                start = regionStart - sequenceOneStartPosition;
                                        }
                                        else
                                        {
                                                start = 0;
                                        }
                                        //right side overlapping exon region
                                        if(regionEnd > sequenceOneEndPosition)
                                        {
                                                end = sequenceOneLength;
                                        }
                                        else
                                        {
                                                end = regionEnd - sequenceOneStartPosition;
                                        }
                                        //cut the region out
                                        string alignmentSeq1Sub = alignmentsWithGapsRemoved[0].substr(start, end);
                                        string alignmentSeq2Sub = alignmentsWithGapsRemoved[1].substr(start, end);
                                        //run calculations on sub strings
                                        compare(alignmentSeq1Sub, alignmentSeq2Sub, IntergenicStats);
                                }
                        }
                        //End of look up for intergenicss *********************************************************************************************************************

			cout << "paircount =  " << pairCount << endl;
			pairCount++;
		}
	}
	inFile.close();

	cout << "**************************************Start of Exons Stats****************************************************************"<< endl;
	printMatches(ExonsStats);
	printGaps(ExonsStats);
	cout << "**************************************End of Exons Stats****************************************************************"<< endl;
	cout << "**************************************Start of Introns Stats****************************************************************"<< endl;
        printMatches(IntronsStats);
        printGaps(IntronsStats);
        cout << "**************************************End of Introns Stats****************************************************************"<< endl;
	cout << "**************************************Start of Codings Stats****************************************************************"<< endl;
        printMatches(CodingsStats);
        printGaps(CodingsStats);
        cout << "**************************************End of Codings Stats****************************************************************"<< endl;
	cout << "**************************************Start of Utr3 Stats****************************************************************"<< endl;
        printMatches(Utr3Stats);
        printGaps(Utr3Stats);
        cout << "**************************************End of Utr3 Stats****************************************************************"<< endl;
	cout << "**************************************Start of Utr5 Stats****************************************************************"<< endl;
        printMatches(Utr5Stats);
        printGaps(Utr5Stats);
        cout << "**************************************End of Utr5 Stats****************************************************************"<< endl;
	cout << "**************************************Start of Promotor length 100 Stats****************************************************************"<< endl;
        printMatches(Promoters100Stats);
        printGaps(Promoters100Stats);
        cout << "**************************************End of Promotor length 100 Stats****************************************************************"<< endl;
	cout << "**************************************Start of Promotor length 200 Stats****************************************************************"<< endl;
        printMatches(Promoters200Stats);
        printGaps(Promoters200Stats);
        cout << "**************************************End of Promotor length 200 Stats****************************************************************"<< endl;
	cout << "**************************************Start of Promotor length 400 Stats****************************************************************"<< endl;
        printMatches(Promoters400Stats);
        printGaps(Promoters400Stats);
        cout << "**************************************End of Promotor length 400 Stats****************************************************************"<< endl;
	cout << "**************************************Start of Promotor length 800 Stats****************************************************************"<< endl;
        printMatches(Promoters800Stats);
        printGaps(Promoters800Stats);
        cout << "**************************************End of Promotor length 800 Stats****************************************************************"<< endl;
	cout << "**************************************Start of Intergenic Stats****************************************************************"<< endl;
        printMatches(IntergenicStats);
        printGaps(IntergenicStats);
        cout << "**************************************End of Intergenic Stats****************************************************************"<< endl;

	
	return 1;
}
