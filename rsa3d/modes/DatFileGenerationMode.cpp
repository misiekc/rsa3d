//
// Created by pkua on 05.09.2019.
//

#include <fstream>

#include "DatFileGenerationMode.h"
#include "../PackingGenerator.h"

DatFileGenerationMode::DatFileGenerationMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("<directory>"));

    this->dirName = positionalArguments[0];
}

double DatFileGenerationMode::getMedian(std::vector<double> &data){
	size_t size = data.size();
	sort(data.begin(), data.end());
	if (size % 2 == 0)
		return (data[size / 2 - 1] + data[size / 2]) / 2;
	else
		return data[size / 2];
}

void DatFileGenerationMode::run() {

    std::string sFile = this->params.getPackingSignature() + ".dat";
    std::ofstream dataFile(sFile);
    if (!dataFile)
    	die("Cannot open file " + sFile + " to store packing info");
    dataFile.precision(std::numeric_limits<double>::digits10 + 1);

    auto filenames = PackingGenerator::findPackingsInDir(dirName);
    std::vector<double> vTimes;
    double sno = 0.0, sno2 = 0.0;
    for (const auto &filename : filenames) {
    	int no1 = lastIndexOf(filename, '_');
    	int no2 = lastIndexOf(filename, '.');
    	std::string seed = filename.substr(no1 + 1, no2 - no1 - 1);

    	Packing packing;
    	packing.restore(filename);

    	dataFile << seed << "\t" << packing.back()->no << "\t" << packing.back()->time << std::endl;
    	std::cout << ".";
    	std::cout.flush();
    	vTimes.push_back(packing.back()->time);
    	sno += (double)packing.back()->no;
    	sno2 += ((double)packing.back()->no)*((double)packing.back()->no);
    }
    std::cout << std::endl;
    double median = getMedian(vTimes);
    /*
        for(size_t i=0; i<vTimes.size(); i++){
        	vTimes[i] = std::abs(vTimes[i]-median);
        }
        double mad = getMedian(vTimes);
    */
    double mad = (vTimes[vTimes.size()/2 + 1] - vTimes[vTimes.size()/2 - 1])/2.0;
    sno /= vTimes.size();
    sno2 /= vTimes.size();

    std::cout << sno << "\t" << std::sqrt((sno2 - sno*sno)/vTimes.size()) << "\t" << median << "\t" << mad << std::endl;
}
