#include <fstream>
#include <vector>
#include <sstream>
#include "mytypes.h"
#include <iomanip>

using namespace MyTypes;

enum class SolverTypes : int {
    WOS, MOI
};
enum class SamplerTypes : int {
    UNIFORM, HALTON
};

struct Config{
    Config(std::string configPath){
        std::ifstream ifile (configPath);
        std::string solverName;
        std::string samplerName;
        ifile >> homeDir >> refName >> expName >> boundName;
        ifile >> numRow >> numCol >> numZ;
        ifile >> stopAt;
        ifile >> samplerName;
        ifile >> solverName;
        ifile >> maxError;
        numCells = numRow * numCol * numZ;
        expName = solverName + "_" + expName;
        if(samplerName == "uniform") samplerType = SamplerTypes::UNIFORM;
        else if(samplerName == "halton") samplerType = SamplerTypes::HALTON;
        else throw std::exception("Unrecognized sampler name");

        if(solverName == "wos") {
            solverType = SolverTypes::WOS;
            if(samplerType == SamplerTypes::HALTON){
                std::cerr << "WoS does not support QMC at the moment, falling back to uniform sampling...";
                samplerType = SamplerTypes::UNIFORM;
                samplerName = "uniform";
            }
        }
        else if(solverName == "moi") {
            solverType = SolverTypes::MOI;
            ifile >> moiNumBounces >> moiNumSamplesPerRay;
        }
        else throw std::exception("Unrecognized solver name");
        std:: ostringstream ss;
        ss << maxError;
        ifile.close();
        info = to_string(numRow) + "x" + to_string(numCol);
        info +="\nBoundary:" + boundName;
        info += "\nEpsilon: " + ss.str();
        info += "\nSampler: " + samplerName;
        info += "\nSolver: " + solverName;
        if(solverType == SolverTypes::MOI){
            info += "\nBounces: " + to_string(moiNumBounces);
            info += "\nSamples/ray: " + to_string(moiNumSamplesPerRay);
        }
    }

    void getBoundary(std::vector<arrayd<2>> &singleSided, std::vector<arrayd<2>> &doubleSided){
        std::ifstream ifile (homeDir+boundName+".txt");
        int numPoints = 0;
        ifile >> numPoints;
        double x,y;
        for(int i = 0; i < numPoints; i++){
            ifile >> x >> y;
            singleSided.emplace_back(x,y);
        }
        ifile >> numPoints;
        for(int i = 0; i < numPoints; i++){
            ifile >> x >> y;
            doubleSided.emplace_back(x,y);
        }
    }

    void getSol(std::unordered_map<int,double> ref) {
        std::ifstream ifile(homeDir + refName + ".txt");
        int numPoints = 0;
        ifile >> numPoints;
        while(numPoints--) {
            int i;
            ifile >> i;
            ifile >> ref[i];
        }
    }
    void saveSol(std::unordered_map<int, double> ref) {
        std::ofstream ofile(homeDir + expName + ".txt");
        for(auto it : ref){
            ofile << it.first << " " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << it.second;
            ofile << "\n";
        }
    }

    void getGrad(arrayd<3> *ref){
        std::ifstream ifile (homeDir+refName+"_grad.txt");
        double tmp;
        double x,y,z;
        for(int i = 0; i < numCells; i++){
            ifile >> x >> y >> z;
            ref[i] = {{x,y,z}};
        }
    }
    void saveGrad(std::vector<arrayd<3>> &ref, int numSample, double error){
        std::ofstream ofile (homeDir+expName+"_grad.txt");
        for(auto x : ref){
            ofile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << x.x();
            ofile << " ";
            ofile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << x.y();
            ofile << " ";
            ofile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << x.z();
            ofile << "\n";
        }
    }
    double maxError;
    int numRow, numCol, numZ, numCells;
    int moiNumBounces, moiNumSamplesPerRay;
    int stopAt;
    std::string expName;
    std::string refName;
    SamplerTypes samplerType;
    SolverTypes solverType;
    std::string homeDir;
    std::string boundName;
    std::string info;
};