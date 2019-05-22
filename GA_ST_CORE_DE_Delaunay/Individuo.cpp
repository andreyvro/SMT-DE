#include "Individuo.h"

Individuo::Individuo() {
    this->grafo = Grafo::getInstancia();
    unsigned int qtdVertices = this->grafo->getQtdVertices();
	// Inicia cromossomo (matriz de adjacência) zerado
    this->matrizAdj.resize(qtdVertices);
    for(unsigned int i = 0; i < qtdVertices; i++) {
        this->matrizAdj.at(i).resize(qtdVertices, false);
    }
}

const unsigned int Individuo::getQtdPtFixo() {
	return this->grafo->getQtdVertices();
}

const unsigned int Individuo::getQtdPtSteiner() {
	return this->qtdPontosSteiner;
}

const unsigned int Individuo::getQtdPtTotal() {
	return this->getQtdPtFixo() + this->qtdPontosSteiner;
}

const unsigned int Individuo::getQtdMaxSteiner() {
	return this->grafo->getQtdVertices() - 2;
}

void Individuo::addPontoSteiner(Vertice vertice) {
	// Adiciona ponto steiner
	this->pontoSteiner.push_back(vertice);
	this->qtdPontosSteiner++;
	// Adiciona no cromossomo
	unsigned int tamCrom = this->matrizAdj.size() + 1;		// Novo tamanho
	this->matrizAdj.resize(tamCrom);
	for(unsigned int i = 0; i < tamCrom; i++) {
        this->matrizAdj[i].resize(tamCrom);
    }
}

void Individuo::remPontoSteiner(unsigned int indice) {
	if (indice >= grafo->getQtdVertices()) {
        // Remove do cromossomo
        this->matrizAdj.erase(this->matrizAdj.begin() + indice);			// Deleta linha
        for(unsigned int i = 0; i < this->matrizAdj.size(); i++) {
            this->matrizAdj[i].erase(this->matrizAdj[i].begin() + indice);	// Deleta coluna
        }
        // Remove ponto steiner
        this->pontoSteiner.erase(this->pontoSteiner.begin() + indice);
        this->qtdPontosSteiner--;
	}
}

Vertice &Individuo::getVertice(unsigned int indice) {
    unsigned int qtdVertices = this->grafo->getQtdVertices();
    if (indice < qtdVertices) {
        return this->grafo->getVertice(indice);
    } else {
        return this->pontoSteiner.at(indice - qtdVertices);
    }
}

double Individuo::getDistancia(unsigned int indice1, unsigned int indice2) {
    unsigned int qtdVerticesGrafo = this->grafo->getQtdVertices();
    if ((indice1 < qtdVerticesGrafo) && (indice2 < qtdVerticesGrafo)) {
        return this->grafo->getDistancia(indice1, indice2);
    } else {
        Vertice &v1 = this->getVertice(indice1);
        Vertice &v2 = this->getVertice(indice2);
        return Aux::getDistancia(v1.getX(), v1.getY(), v2.getX(), v2.getY());
    }
}

double Individuo::getAngulo(unsigned int indice1, unsigned int indice2, unsigned int indice3) {
    Vertice &v1 = this->getVertice(indice1);
    Vertice &v2 = this->getVertice(indice2);
    Vertice &v3 = this->getVertice(indice3);
    return Aux::getAngulo(v1.getX(), v1.getY(), v2.getX(), v2.getY(), v3.getX(), v3.getY());
}

void Individuo::calcularFitness() {
    // Soma distância euclidiana dos vertices que posuem ligação
	double fitness = 0;
	unsigned int totalVertices = this->getQtdPtTotal();
	#pragma omp parallel for
	for (unsigned int i = 0; i < totalVertices; i++) {			// Linha da meia matriz de adjacência
		#pragma omp parallel for reduction(+:fitness)
		for (unsigned int j = i + 1; j < totalVertices; j++) {	// Coluna da meia matriz de adjancência
			if (this->getGene(i, j) == true) {		            // Se o gene for true, os vertices estão ligados
				fitness += this->getDistancia(i, j);            // Calcula distancia euclidiana entre dois vértices
            }
		}
	}
	this->setFitness(fitness);                                  // Aplica valor de fitness ao indivíduo
}

void Individuo::setGene(unsigned int x, unsigned int y, bool valor) {
    // Seta valor no cromossomo
	this->matrizAdj[x][y] = valor;
	this->matrizAdj[y][x] = valor;
}

const bool Individuo::getGene(unsigned int x, unsigned int y) {
	return this->matrizAdj[x][y];
}

unsigned int Individuo::getVertGrau(unsigned int vertice) {
    unsigned int tamCrom = this->matrizAdj.size();
    unsigned int grau = 0;
    for (unsigned int col = 0; col < tamCrom; col++) {
        if (this->matrizAdj[vertice][col] == 1) {
            grau++;
        }
    }
    return grau;
}

void Individuo::getPtsAdjacentes(unsigned int vertice, vector<unsigned int> &adj) {
	unsigned int tamCrom = this->matrizAdj.size();
    for (unsigned int col = 0; col < tamCrom; col++) {
        if (this->matrizAdj[vertice][col] == true) {
            adj.push_back(col);
        }
    }
}

void Individuo::setFitness(double fitness) {
	this->fitness = fitness;
}

const double Individuo::getFitness() {
	return this->fitness;
}

const double Individuo::getXMin() {
    return this->grafo->getXMin();
}
const double Individuo::getXMax() {
    return this->grafo->getXMax();
}
const double Individuo::getYMin() {
    return this->grafo->getYMin();
}
const double Individuo::getYMax() {
    return this->grafo->getYMax();
}
