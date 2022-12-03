#include <iostream>
#include <cmath>
#include <vector>
#include "pbplots.hpp"
#include "supportLib.hpp"
using namespace std;

class F1Model {
private:
    // Coefficients
    float tireDegVD;
    float fuelBurnG;
    float avgVelV;
    float fuelBurnV;
    float tireDegP;
    float tireDegDV;
    float tireDegDP;
    float tireTempP;
    float avgVelP;
    float tireTempD;
    float tireGripPG;
    float lapTimeL;
    float laptimeLG;
    float degTempT;
    float tireDegVP;
    // Attributes
    float tireDegradation = 0.;
    float tireTemperature = 75.;
    float tireGrip = 1.;
    float averageVelocity = 10; // has to start at non zero because there's no description of getting up to speed, consider an acc param
    float amountOfFuel = 100;

    float stepSize = 0.1;


public:
    F1Model(float a, float b, float c, float d, float e, float f, float g, float h, float i, float j, float k, float l, float m, float n, float o) {
        tireDegVD = a;
        fuelBurnG = b;
        avgVelV = c;
        fuelBurnV = d;
        tireDegP = e;
        degTempT = f;
        tireDegDV = g;
        tireDegDP = h;
        tireTempP = i;
        avgVelP = j;
        tireTempD = k;
        tireGripPG = l;
        lapTimeL = m;
        laptimeLG = n;
        tireDegVP = o;

    }
    float simulateIter();
    // Differential functions arranged by hierarchy
    float differentialAverageVelocity(float V, float T, float alpha); // 1
    float differentalVelocityRungeKutta(float V, float T, float alpha); //1 RK
    float differentialTireDegradation(float T, float V, float P); // 2
    float differentalTireDegradationRungeKutta(float T, float V, float P); // 2 RK
    float differentialTireTemperature(float T); // 3
    float differentialTireTemperatureRungeKutta(float T); // 3 RK
    float lapTime(float distance, float velocity);
};


int main(){
    F1Model f1Model(0.01, 0.01, 0.1, 0.01, 0.01, 0.01, 0.0001, 0.0001, 0.01, 0.01, 0.05, 0.01, 0.01, 0.01, 0.0001);
    RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();
    StringReference *errorMessage = CreateStringReferenceLengthValue(0, L' ');

    bool success;

    float laptime;
    vector<double> timestep;
    vector<double> lapTimeVec = {};
    for(int t  = 0; t < 100; t++){
        laptime = f1Model.simulateIter();
        lapTimeVec.push_back(laptime);
        timestep.push_back(t);
        cout << "laptime: " << laptime << endl;
    }

    success = DrawScatterPlot(imageReference, 600, 400, &timestep, &lapTimeVec, errorMessage);

    if(success){
        vector<double> *pngdata = ConvertToPNG(imageReference->image);
        WriteToFile(pngdata, "example1.png");
        DeleteImage(imageReference->image);
    }else{
        cerr << "Error: ";
        for(int i = 0; i < errorMessage->string->size(); i++){
            wcerr << errorMessage->string->at(i);
        }
        cerr << endl;
    }
    FreeAllocations();
}

float F1Model::simulateIter() {
    // Calculate differentials
    float deltaVel = differentalVelocityRungeKutta(averageVelocity, tireDegradation, amountOfFuel);
    float deltaDeg = differentalTireDegradationRungeKutta(tireDegradation, averageVelocity, tireTemperature);
    float deltaTemp = differentialTireTemperatureRungeKutta(tireDegradation);
    // Update attributes
    tireDegradation += deltaDeg;
    tireTemperature += deltaTemp;
    averageVelocity += deltaVel;

    cout << "tire deg: " << tireDegradation << endl;
    cout << "tire temp: " << tireTemperature << endl;
    cout << "avg vel: " << averageVelocity << endl << endl;

    // Calculate lap time
    float distance = 100;
    float lapTimeAvg = lapTime(distance, averageVelocity);

    return lapTimeAvg;
}


float F1Model::differentialAverageVelocity(float V, float T, float alpha) {

    return avgVelV * V - tireDegVD * T * averageVelocity + fuelBurnV / alpha;
}

float F1Model::differentalVelocityRungeKutta(float V, float T, float alpha){

    float k1 = differentialAverageVelocity(V, T, alpha); //N term
    float k2 = differentialAverageVelocity(V + stepSize*k1/2, T, alpha);
    float k3 = differentialAverageVelocity(V + stepSize*k2/2, T, alpha);
    float k4 = differentialAverageVelocity(V + stepSize*k3, T, alpha);

    return stepSize*(k1 + 2*k2 + 2*k3 + k4)/6;

}

float F1Model::differentialTireDegradation(float T, float V, float P) {

    return tireDegDV * T * V + tireDegDP * T * P + V * P * tireDegVP;
}

float F1Model::differentalTireDegradationRungeKutta(float T, float V, float P){

    float k1 = differentialTireDegradation(T, V, P); //N term
    float k2 = differentialTireDegradation(T + stepSize*k1/2, V, P);
    float k3 = differentialTireDegradation(T + stepSize*k2/2, V, P);
    float k4 = differentialTireDegradation(T + stepSize*k3, V, P);

    return stepSize*(k1 + 2*k2 + 2*k3 + k4)/6;
}

float F1Model::differentialTireTemperature(float T) {

    return tireTempD * sqrt(T);
}

float F1Model::differentialTireTemperatureRungeKutta(float T){

    float k1 = differentialTireTemperature(T); //N term
    float k2 = differentialTireTemperature(T + stepSize*k1/2);
    float k3 = differentialTireTemperature(T + stepSize*k2/2);
    float k4 = differentialTireTemperature(T + stepSize*k3);

    return stepSize*(k1 + 2*k2 + 2*k3 + k4)/6;

}

float F1Model::lapTime(float distance, float velocity){
    float lapTime = distance / velocity;
    return lapTime;
}

