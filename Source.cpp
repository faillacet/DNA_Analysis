#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <time.h>
using namespace std;

class Calculations {
private:
	string *cheese;					//dynamicly allocated arrays
	int *sums;
	double *SD;
	char *c;
	int size = 0;					//results storage
	double mean = 0;
	double standardDeviation = 0;
	double variance = 0;
	int totalSum = 0;
	int countA = 0;					//count of each nucleotide
	int countC = 0;
	int countT = 0;
	int countG = 0;
	double probA = 0;				//frequency of each nucleotide
	double probC = 0;
	double probT = 0;
	double probG = 0;
	int bigram[16];					//number of times each bigram occurs
	float bigramFrequency[16];		//frequency of each bigram
public:
	void getFile(string file); //function interates through all strings in sample and gets length of each string and stores them in an array
	void copyFile(string file);//used in get file to actually copy the data to an array
	void calcAll();				//all calculations done in this function
	int calcSum(string x);		//calc sum of DNA strings
	double calcMean();			//calc mean of DNA strings
	double calcVar();			//calc variance of DNA strings
	double calcSD();			//calc standard deviation of length of DNA string
	void calcProbability();		//calc relative probabilities of each nucleotide
	void bigramScan(string x);	//calc number of times each bigram occurs
	void calcBigramFrequency();	//calc relative probabilities of bigrams
	void calcTotalSum();

	double gaussian();			//uses some math stuff to create a normal distruibution
	string generateString();	//generate string based on normal distribution (of calcd data)
	void outputFile();

	void runProgram();			//USES everything else to get user input for file path then preforms all operations needed
	bool runAgain();			//returns true as long as user wants to do another to rerun the program, if false program exits

	Calculations() {
		srand((unsigned)time(NULL));		//code from: https://www.softwaretestinghelp.com/random-number-generator-cpp/#:~:text=We%20can%20use%20srand%20(),value%20either%20float%20or%20double.
	};										//used to seed the RNG

	~Calculations() {
		delete[] cheese;
		delete[] sums;
		delete[] SD;
		cout << "Dynamically allocated memory wiped." << endl;
		cout << "See you later!" <<endl;
	}
};

double Calculations::gaussian() {		//Calculations for the distribution
	double pi = 2 * acos(0.0);		//idea to calc pi from geeksforgeeks
	double a, b, c, d;
	a = (double)rand() / RAND_MAX;	//idea from ^^
	b = (double)rand() / RAND_MAX;	//^^			-- not exactly random cuz it is seeded by time but good enough for our purposes
	c = sqrt(-2 * log(a))*cos(2 * b*pi);
	d = standardDeviation * c + mean;
	return d;
}

string Calculations::generateString() {
	int stringLength = round(gaussian());
	c = new char[stringLength];
	double a;
	for (int i = 0; i < stringLength; ++i) {
		a = (double)rand() / RAND_MAX;
		if (a < probA)
			c[i] = 'A';
		else if (a < probA + probC)
			c[i] = 'C';
		else if (a < probA + probC + probT)
			c[i] = 'T';
		else if (a <= probA + probC + probG + probT)
			c[i] = 'G';
	}
	string s = "";
	for (int i = 0; i < stringLength; ++i)
		s = s + c[i];
	delete[] c;
	return s;
}

void Calculations::calcTotalSum() {
	for (int i = 0; i < size; ++i) {
		totalSum += sums[i];
	}
}

void Calculations::getFile(string file) {
	ifstream myFile(file);
	if (myFile.is_open()) {
		int count = 0;
		string x;
		while (getline(myFile, x)) {
			count++;
		}
		size = count;
		cheese = new string[count];
		sums = new int[count];
		SD = new double[count];
		copyFile(file);
	}
	else {
		cout << "Error opening file, exiting program..." << endl;
	}
}

void Calculations::copyFile(string file) {
	ifstream myFile(file);
	if (myFile.is_open()) {
		int count = 0;
		string x;
		while (getline(myFile, x)) {
			cheese[count] = x;
			count++;
		}
		myFile.close();
	}
	else
		cout << "Error copying file, exiting program..." << endl;
}

int Calculations::calcSum(string x) {
	int sum = 0;
	for (int i = 0; i < x.length(); ++i) {
		if (x[i] == 'A' || x[i] == 'a') {
			sum += 1;
			countA += 1;
		}
		else if (x[i] == 'C' || x[i] == 'c') {
			sum += 1;
			countC += 1;
		}
		else if (x[i] == 'T' || x[i] == 't') {
			sum += 1;
			countT += 1;
		}
		else if (x[i] == 'G' || x[i] == 'g') {
			sum += 1;
			countG += 1;
		}
	}
	return sum;
}

double Calculations::calcMean() {
	double map = 0;
	for (int i = 0; i < size; ++i) {
		map += sums[i];
	}
	map = map / size;
	return map;
}

double Calculations::calcSD() {
	return sqrt(variance);
}

double Calculations::calcVar() {
	for (int i = 0; i < size; ++i) {
		SD[i] = (sums[i] - mean)*(sums[i] - mean);
	}
	double map = 0;
	for (int i = 0; i < size; ++i) {
		map += SD[i];
	}
	map = map / (size - 1);
	return map;
}

void Calculations::calcProbability() {
	double total = countA + countC + countT + countG;
	probA = countA / total;
	probC = countC / total;
	probT = countT / total;
	probG = countG / total;
}

void Calculations::calcBigramFrequency() {
	float total = 0;
	for (int i = 0; i < 16; ++i) {
		total += bigram[i];
	}
	for (int i = 0; i < 16; ++i) {
		bigramFrequency[i] = bigram[i] / total;
	}
}

void Calculations::bigramScan(string x) {													//is there an easier way to do this?
	for (int i = 0; i < x.length(); i += 2) {
		int z = i + 1;
		if ((x[i] == 'A' || x[i] == 'a') && (x[i + 1] == 'A' || x[i + 1] == 'a')) {
			bigram[0] += 1;
		}
		if ((x[i] == 'A' || x[i] == 'a') && (x[i + 1] == 'C' || x[i + 1] == 'c'))
			bigram[1] += 1;
		if ((x[i] == 'A' || x[i] == 'a') && (x[i + 1] == 'T' || x[i + 1] == 't'))
			bigram[2] += 1;
		if ((x[i] == 'A' || x[i] == 'a') && (x[i + 1] == 'G' || x[i + 1] == 'g'))
			bigram[3] += 1;
		if ((x[i] == 'C' || x[i] == 'c') && (x[i + 1] == 'A' || x[i + 1] == 'a'))
			bigram[4] += 1;
		if ((x[i] == 'C' || x[i] == 'c') && (x[i + 1] == 'C' || x[i + 1] == 'c'))
			bigram[5] += 1;
		if ((x[i] == 'C' || x[i] == 'c') && (x[i + 1] == 'T' || x[i + 1] == 't'))
			bigram[6] += 1;
		if ((x[i] == 'C' || x[i] == 'c') && (x[i + 1] == 'G' || x[i + 1] == 'g'))
			bigram[7] += 1;
		if ((x[i] == 'T' || x[i] == 't') && (x[i + 1] == 'A' || x[i + 1] == 'a'))
			bigram[8] += 1;
		if ((x[i] == 'T' || x[i] == 't') && (x[i + 1] == 'C' || x[i + 1] == 'c'))
			bigram[9] += 1;
		if ((x[i] == 'T' || x[i] == 't') && (x[i + 1] == 'T' || x[i + 1] == 't'))
			bigram[10] += 1;
		if ((x[i] == 'T' || x[i] == 't') && (x[i + 1] == 'G' || x[i + 1] == 'G'))
			bigram[11] += 1;
		if ((x[i] == 'G' || x[i] == 'g') && (x[i + 1] == 'A' || x[i + 1] == 'a'))
			bigram[12] += 1;
		if ((x[i] == 'G' || x[i] == 'g') && (x[i + 1] == 'C' || x[i + 1] == 'c'))
			bigram[13] += 1;
		if ((x[i] == 'G' || x[i] == 'g') && (x[i + 1] == 'T' || x[i + 1] == 't'))
			bigram[14] += 1;
		if ((x[i] == 'G' || x[i] == 'g') && (x[i + 1] == 'G' || x[i + 1] == 'g'))
			bigram[15] += 1;
	}
	int y = x.length() % 2;
	int z = x.length() - 1;
	if (y == 1) {																			//if odd compare first to last :(
		if ((x[0] == 'A' || x[0] == 'a') && (x[z] == 'A' || x[z] == 'a'))
			bigram[0] += 1;
		else if ((x[0] == 'A' || x[0] == 'a') && (x[z] == 'C' || x[z] == 'c'))
			bigram[1] += 1;
		else if ((x[0] == 'A' || x[0] == 'a') && (x[z] == 'T' || x[z] == 't'))
			bigram[2] += 1;
		else if ((x[0] == 'A' || x[0] == 'a') && (x[z] == 'G' || x[z] == 'g'))
			bigram[3] += 1;
		else if ((x[0] == 'C' || x[0] == 'c') && (x[z] == 'A' || x[z] == 'a'))
			bigram[4] += 1;
		else if ((x[0] == 'C' || x[0] == 'c') && (x[z] == 'C' || x[z] == 'c'))
			bigram[5] += 1;
		else if ((x[0] == 'C' || x[0] == 'c') && (x[z] == 'T' || x[z] == 't'))
			bigram[6] += 1;
		else if ((x[0] == 'C' || x[0] == 'c') && (x[z] == 'G' || x[z] == 'g'))
			bigram[7] += 1;
		else if ((x[0] == 'T' || x[0] == 't') && (x[z] == 'A' || x[z] == 'a'))
			bigram[8] += 1;
		else if ((x[0] == 'T' || x[0] == 't') && (x[z] == 'C' || x[z] == 'c'))
			bigram[9] += 1;
		else if ((x[0] == 'T' || x[0] == 't') && (x[z] == 'T' || x[z] == 't'))
			bigram[10] += 1;
		else if ((x[0] == 'T' || x[0] == 't') && (x[z] == 'G' || x[z] == 'G'))
			bigram[11] += 1;
		else if ((x[0] == 'G' || x[0] == 'g') && (x[z] == 'A' || x[z] == 'a'))
			bigram[12] += 1;
		else if ((x[0] == 'G' || x[0] == 'g') && (x[z] == 'C' || x[z] == 'c'))
			bigram[13] += 1;
		else if ((x[0] == 'G' || x[0] == 'g') && (x[z] == 'T' || x[z] == 't'))
			bigram[14] += 1;
		else if ((x[0] == 'G' || x[0] == 'g') && (x[z] == 'G' || x[z] == 'g'))
			bigram[15] += 1;
	}
}

void Calculations::outputFile() {
	ofstream myFile("TrentonFaillace.txt");
	if (myFile.is_open()) {
		myFile << "Trenton Faillace" << endl;
		myFile << "ID: 2382418" << endl << endl;
		myFile << "Statistics Summary: " << endl;
		myFile << "Total Sum of all strings: " << totalSum << endl;
		myFile << "Mean: " << mean << endl;
		myFile << "Standard Deviation: " << standardDeviation << endl;
		myFile << "Variance: " << variance << endl << endl;
		myFile << "Relative frequency of nucleotides in decimal:" << endl;
		myFile << "A: " << probA << endl;
		myFile << "C: " << probC << endl;
		myFile << "T: " << probT << endl;
		myFile << "G: " << probG << endl << endl;
		myFile << "Relative frequency of bigrams in decimal:" << endl;
		myFile << "AA: " << bigramFrequency[0] << endl;
		myFile << "AC: " << bigramFrequency[1] << endl;
		myFile << "AT: " << bigramFrequency[2] << endl;
		myFile << "AG: " << bigramFrequency[3] << endl;
		myFile << "CA: " << bigramFrequency[4] << endl;
		myFile << "CC: " << bigramFrequency[5] << endl;
		myFile << "CT: " << bigramFrequency[6] << endl;
		myFile << "CG: " << bigramFrequency[7] << endl;
		myFile << "TA: " << bigramFrequency[8] << endl;
		myFile << "TC: " << bigramFrequency[9] << endl;
		myFile << "TT: " << bigramFrequency[10] << endl;
		myFile << "TG: " << bigramFrequency[11] << endl;
		myFile << "GA: " << bigramFrequency[12] << endl;
		myFile << "GC: " << bigramFrequency[13] << endl;
		myFile << "GT: " << bigramFrequency[14] << endl;
		myFile << "GG: " << bigramFrequency[15] << endl << endl;
		string x;
		for (int i = 0; i < 1000; ++i) {
			x = generateString();
			myFile << x << endl;
		}
		myFile.close();
	}
	else {
		cout << "Failed to write to file" << endl;
	}
}

void Calculations::calcAll() {
	int map = 0;
	for (int i = 0; i < size; ++i) {
		map = calcSum(cheese[i]);
		sums[i] = map;
	}

	mean = calcMean();
	variance = calcVar();
	standardDeviation = calcSD();
	calcProbability();
	calcTotalSum();

	for (int i = 0; i < 16; ++i) {
		bigram[i] = 0;
	}
	for (int i = 0; i < size; ++i) {
		bigramScan(cheese[i]);
	}
	calcBigramFrequency();
}

void Calculations::runProgram() {
	string x;
	cout << "Please enter name of file you would like to process.  File must be in same directory as this program." << endl;
	cout << "EXAMPLE: DNA.txt" << endl;
	cout << "File Name: ";
	cin >> x;

	getFile(x);

	calcAll();
	outputFile();
}

bool Calculations::runAgain() {
	char userInput;
	cout << "Would you like to process another list? (Type Y/N)" << endl;
	cin >> userInput;
	if (userInput == 'N' || userInput == 'n') {
		return false;
	}
	else if (userInput == 'Y' || userInput == 'y') {
		return true;
	}
	else {
		cout << "Error, response is incorrect." << endl;
		cout << "Exiting program..." << endl;
		return false;
	}
}

int main(int argc, char* argv[]) {								//main is super small plz give me good grade
	Calculations calc;
	bool flag = true;
	while (flag == true) {
		calc.runProgram();
		flag = calc.runAgain();
	}
	return 0;
}
