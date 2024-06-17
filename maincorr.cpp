#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <iomanip>

// Función calcular vecinos
inline int mod(int a, int b) {
    return (a % b + b) % b;
}

// Clase para el modelo XY
class XYModel {
public:
    XYModel(int L, double temp, int seed);
    void initialize();
    double wolffClusterUpdate(int niters);
    void saveConfiguration(const std::string &filename);
    double calculateEnergy() const;
    void saveEnergy(const std::string &filename, double energy);
    std::vector <double> calculateCorr();
    void saveCorr(const std::string &filename, std::vector<double> &Corr);
    std::mt19937 gen;
    double temp;
    std::vector<double> spins;

private:
    int L;
    int N;
    std::vector<int> neighbors;
    std::uniform_real_distribution<> dis;
};

// Constructor
XYModel::XYModel(int L, double temp, int seed)
    : L(L), temp(temp), N(L * L), gen(seed), dis(0.0, 2 * M_PI) {
    initialize();
}

// Inicialización de la configuración y vecinos
void XYModel::initialize() {
    spins.resize(N);
    neighbors.resize(4 * N);
    for (int i = 0; i < N; ++i) {
        spins[i] = dis(gen);
        int x = i / L;
        int y = i % L;
        neighbors[4 * i + 0] = mod(x - 1, L) * L + y; // izquierda
        neighbors[4 * i + 1] = mod(x + 1, L) * L + y; // derecha
        neighbors[4 * i + 2] = x * L + mod(y - 1, L); // arriba
        neighbors[4 * i + 3] = x * L + mod(y + 1, L); // abajo
    }
}

// Actualización del sistema por algoritmo Wolff cluster
double XYModel::wolffClusterUpdate(int niters) {
    std::uniform_real_distribution<> rand_prob(0.0, 1.0);
    double beta = 1.0 / temp;
    double avg_cluster_size = 0.0;
    //niters = std::ceil(niters / 4) + 1;                      //Descomentar esta línea para mayor eficiencia

    for (int n = 0; n < niters; ++n) {
        std::vector<int> cluster(N, 0);
        std::vector<int> stack(N + 1, 0);

        int site = rand_prob(gen) * N;
        double rand_angle = dis(gen);
        cluster[site] = 1;
        stack[0] = site + 1;

        int sc_in = 0, sc_out = 1;

        while (sc_in < sc_out) {
            int current_site = stack[sc_in++] - 1;
            double spin_angle = spins[current_site];
            spins[current_site] = 2 * rand_angle - spin_angle;
            avg_cluster_size++;

            for (int k = 0; k < 4; ++k) {
                int neighbor = neighbors[4 * current_site + k];
                if (!cluster[neighbor]) {
                    double neighbor_angle = spins[neighbor];
                    double dE = cos(spin_angle - neighbor_angle) - cos(2 * rand_angle - spin_angle - neighbor_angle);
                    if (rand_prob(gen) < 1.0 - exp(-beta * dE)) {
                        cluster[neighbor] = 1;
                        stack[sc_out++] = neighbor + 1;
                    }
                }
            }
        }
    }

    return avg_cluster_size / (niters);
}

// Calcular la energía de la configuración
double XYModel::calculateEnergy() const {
    double energy = 0.0;
    for (int i = 0; i < N; ++i) {
        double spin_angle = spins[i];
        for (int k = 0; k < 4; ++k) {
            int neighbor = neighbors[4 * i + k];
            double neighbor_angle = spins[neighbor];
            energy += -cos(spin_angle - neighbor_angle);
        }
    }
    return energy;
}

// Guardar energía en un archivo
void XYModel::saveEnergy(const std::string &filename, double energy) {
    std::ofstream file(filename, std::ios::app);
    file << energy << "\n";
    file.close();
}



// Calcular la función de correlación
std::vector<double> XYModel::calculateCorr() {
	std::uniform_real_distribution<> rand_prob(0.0, 1.0);
	std::vector<double> Corr (L-1,0);
	int file = rand_prob(gen)*(L-1);
	for ( int site = 0; site < L ; ++site){
		double current_angle = spins[ mod(site , L)*L + file ];
		for ( int neigh = 1; neigh < L; ++neigh ){
			int neighbor_site = mod(site + neigh, L)*L + file;
			double neighbor_angle = spins[neighbor_site];
			Corr[neigh-1] += cos(current_angle - neighbor_angle) / L;	
		}
	}
	return Corr;
}	


// Guardar correlaciones en un archivo
void XYModel::saveCorr(const std::string &filename, std::vector<double> &Corr) {
    std::ofstream file(filename, std::ios::app); 
    for (const auto &corr : Corr) {
        file << corr << " ";
    }
    file << "\n"; 
    file.close();
}
	
	

// Guardar configuración en un archivo
void XYModel::saveConfiguration(const std::string &filename) {
    std::ofstream file(filename); 
    for (const auto &spin : spins) {
        file << spin << " ";
    }
    file << "\n"; 
    file.close();
}

// Estructura para las réplicas
struct Replica {
    XYModel model;
    double temperature;
    double energy;
    int niters;
    Replica(int L, double temp, int seed) : model(L, temp, seed), temperature(temp), energy(0), niters(0) {}
};

// Función principal
int main(int argc, char* argv[]) {
    if (argc != 9) {
        std::cerr << "Usage: " << argv[0] << " <L> <Tmin> <Tmax> <steps> <output_dir> <therm> <meas> <mod>" << std::endl;  //se asegura de que se  han introducido correctamente los inputs
        return 1;
    }

    int L = std::stoi(argv[1]);            //Longitud de la red
    double Tmin = std::stod(argv[2]);      //Rango de temperaturas
    double Tmax = std::stod(argv[3]);
    int steps = std::stoi(argv[4]);
    std::string output_dir = argv[5];      //Nombre del directorio en el que se guardarán los outputs
    int therm = std::stoi(argv[6]);        //Bins de termalización
    int meas = std::stoi(argv[7]);         //Bins de medidción
    int mod = std::stoi(argv[8]);          //Parámetro para modo medida o no medida. Si == 1 se calcularán y guardarán las funciones de correlación 
    
    int PTparam = 0;                       //Parámetro para Parallel Tempering. Sirve para intercalar las propuestas entre temperaturas pares e impares
    

    std::vector<double> temperatures;
    double dT = (Tmax - Tmin) / (steps - 1); 
    for (int i = 0; i < steps; ++i) {
        temperatures.push_back(Tmin + i * dT);
    }

    // Crear directorio de salida
    if (access(output_dir.c_str(), F_OK) == -1) {
        mkdir(output_dir.c_str(), 0777);
    }

    // Guardar temperaturas en variables.data
    std::ofstream variables_file(output_dir + "/variables.data");
    for (const auto &temp : temperatures) {
        variables_file << temp << " ";
    }
    variables_file.close();

    // Inicializar réplicas
    std::vector<Replica> replicas;
    for (const auto& temp : temperatures) {
        replicas.emplace_back(L, temp, 42);
    }

    // Pretermalización para encontrar el número adecuado de iteraciones (niters) para cada réplica
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < replicas.size(); ++i) {
        double avg_cluster_size = 0.0;
        int pretherm_steps = 1000;
        for (int step = 0; step < pretherm_steps; ++step) {
            avg_cluster_size += replicas[i].model.wolffClusterUpdate(1);
            }
        
        avg_cluster_size /= pretherm_steps;
        
        replicas[i].niters = std::ceil((L * L) / avg_cluster_size);
        std::cout << "For T=" << replicas[i].temperature << ":  niters = " << replicas[i].niters << std::endl;
    	
    }
    

    int therm_steps = therm*100; // Número de pasos de termalización
    int measure_steps = meas*100; // Número de pasos de medición

    //Termalización y medición con PT
    for (int step = 0; step < therm_steps + measure_steps; ++step) {
        //Termalización
        if (step < therm_steps){
        
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < replicas.size(); ++i) {
            replicas[i].model.wolffClusterUpdate(replicas[i].niters);
            replicas[i].energy = replicas[i].model.calculateEnergy();
        }
	if ( step < therm_steps and step % 100 == 0 ){
	std::cout << "Intento de PT en termalización" << step/100 << "de" << therm_steps/100 << std::endl;
	PTparam=1-PTparam;
        // Propuesta de intercambio de réplicas
        for (size_t i = PTparam; i < replicas.size() - 1; i += 2) {
            double dE = replicas[i + 1].energy - replicas[i].energy;
            double dBeta = 1.0 / replicas[i + 1].temperature - 1.0 / replicas[i].temperature;
            double acceptance_ratio = std::exp(-dBeta * dE);
            if (std::uniform_real_distribution<>(0.0, 1.0)(replicas[i].model.gen) < acceptance_ratio) {
                std::swap(replicas[i], replicas[i + 1]);
            	std::vector<double> spinsi = replicas[i].model.spins;
            	std::vector<double> spinsj = replicas[i+1].model.spins;
            	replicas[i].model.spins = spinsj;
            	replicas[i+1].model.spins = spinsi;
            }
        }
        }
        }

        // Medición después de la termalización
if (step >= therm_steps) {
    double avg_cluster_size = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:avg_cluster_size)
    for (size_t i = 0; i < replicas.size(); ++i) {
    	avg_cluster_size += replicas[i].model.wolffClusterUpdate(replicas[i].niters);
        replicas[i].energy = replicas[i].model.calculateEnergy(); 
    	if ( mod == 1){
    		// Guardar la energía
        	std::ostringstream energy_file;
        	energy_file << output_dir << "/EnergyatT=" << std::setw(5) << std::setfill('0') << static_cast<int>(replicas[i].model.temp * 10000) << ".data";
        	replicas[i].model.saveEnergy(energy_file.str(), replicas[i].energy);
    	
    		//Guardar la correlación
        	std::vector<double> Corr = replicas[i].model.calculateCorr();
        	std::ostringstream corr_file;
            	corr_file << output_dir << "/corratT=" << std::setw(5) << std::setfill('0') << static_cast<int>(replicas[i].model.temp * 10000) << ".data";
            	replicas[i].model.saveCorr(corr_file.str() , Corr );
        	
    	}
    }
    if ( step % 100 == 0 ){
    std::cout << "Intento de PT en medida" << step/100 << "de" << (therm_steps+measure_steps)/100 << std::endl;
    PTparam=1-PTparam;
    // Propuesta de intercambio de réplicas
    for (size_t i = PTparam; i < replicas.size() - 1; i += 2) {
        double dE = replicas[i + 1].energy - replicas[i].energy;
        double dBeta = 1.0 / replicas[i + 1].temperature - 1.0 / replicas[i].temperature;
        double acceptance_ratio = std::exp(-dBeta * dE);
        if (std::uniform_real_distribution<>(0.0, 1.0)(replicas[i].model.gen) < acceptance_ratio) {
            std::swap(replicas[i], replicas[i + 1]);
            std::vector<double> spinsi = replicas[i].model.spins;
            std::vector<double> spinsj = replicas[i+1].model.spins;
            replicas[i].model.spins = spinsj;
            replicas[i+1].model.spins = spinsi;
            
        }
    }
}
}

    // Guardar la última configuración después de todos los pasos de termalización y medición
    if (step == therm_steps + measure_steps - 1) {
        #pragma omp parallel for
        for (size_t i = 0; i < replicas.size(); ++i) {
            std::ostringstream temp_file;
            temp_file << output_dir << "/configatT=" << std::setw(5) << std::setfill('0') << static_cast<int>(replicas[i].model.temp * 10000) << ".data";
            replicas[i].model.saveConfiguration(temp_file.str());
        }
    }
}
return 0;
}





