#ifndef GAAMS_H
#define GAAMS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include "Populacao.h"
#include "Individuo.h"
#include "Aux.h"

using namespace std;

class GA {
	private:
		Populacao *populacao;                                                   // População
		float indCruzamento;                                                    // Índice de Cruzamento
		float indMutacao;                                                       // Índice de Mutação
		double fitAtingir;
		string log;                                                             // Log de execução
        vector< Vertice > ptSt;                                                 // Possíveis pontos Steiner iniciais

		void criarArvoreKruskal(Individuo &ind);

		void getSelecao(vector<unsigned int> &ind, unsigned int qtd = 3, const unsigned int naoSel = -1);
		void setDE_rand_1(unsigned int iInd, Individuo &filho, float f1, float cr);                                     // DE/rand/1
		void setDE_best_1(unsigned int iInd, Individuo &filho, float f1, float cr, unsigned int best);                  // DE/best/1
		void setDE_rand_2(unsigned int iInd, Individuo &filho, float f1, float f2, float cr);                           // DE/rand/2
		void setDE_best_2(unsigned int iInd, Individuo &filho, float f1, float f2, float cr, unsigned int best);        // DE/best/2
		void setDE_curToRand_1(unsigned int iInd, Individuo &filho, float f1, float f2, float cr);                      // DE/ctor/1
		void setDE_curToBest_2(unsigned int iInd, Individuo &filho, float f1, float f2, float cr, unsigned int best);   // DE/ctob/1

		unsigned int addPontoSteiner(Individuo &ind, unsigned int v1, unsigned int v2, unsigned int v3);            // Insere ponto Steiner através de v1-v2, v1-v3
		void remPontoSteiner(Individuo &ind, unsigned int s1, unsigned int v1, unsigned int v2, unsigned int v3);   // Remove ponto Steiner e une v1-v2, v1-v3
        void criarArvoreSteiner(Individuo &ind);                                                                    // Cria árvore de Steiner através de árvore geradora mínima

		void reposicionaPontosSmith(Individuo &individuo, double prec = 0.001);
		double length(Individuo &ind, vector< vector<unsigned int> > &lstAdj, vector< vector< double > > &lstDist);
		double error(Individuo &ind, vector< vector<unsigned int> > &lstAdj, vector< vector< double > > &lstDist);
		void optimize(Individuo &ind, vector< vector<unsigned int> > &lstAdj, vector< vector< double > > &lstDist, double tol);

struct DisjointSets {
    int *parent, *rnk;
    int n;

    // Constructor.
    DisjointSets(int n) {
        // Allocate memory
        this->n = n;
        parent = new int[n + 1];
        rnk = new int[n + 1];

        // Initially, all vertices are in different sets and have rank 0.
        for (int i = 0; i <= n; i++) {
            rnk[i] = 0;
            parent[i] = i;  //every element is parent of itself
        }
    }

    // Find the parent of a node 'u'
    // Path Compression
    int find(int u) {
        //Make the parent of the nodes in the path from u--> parent[u] point to parent[u]
        if (u != parent[u])
            parent[u] = find(parent[u]);
        return parent[u];
    }

    // Union by rank
    void merge(int x, int y) {
        x = find(x), y = find(y);

        //Make tree with smaller height a subtree of the other tree
        if (rnk[x] > rnk[y])
            parent[y] = x;
        else // If rnk[x] <= rnk[y]
            parent[x] = y;

        if (rnk[x] == rnk[y])
            rnk[y]++;
    }
};

	public:
		GA(Populacao *pop, unsigned int qtdInd, float indCruz, float indMut, double fitAt);   // Cria uma população inicial
		void iniciarGA(unsigned int numIteracoes, bool exibirInfo);                                                 // Inicia Algoritmo
		string getLog();                                                                                            // Retorna melhor fitness das gerações
};

#endif // GAAMS_H
