#ifndef GMSTClustering_h
#define GMSTClustering_h 1

#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <vector>
#include "../TGlauber/TGlauberMC.hh"
#include "../TGlauber/TGlauNucleon.hh"
#include "TObject.h"
#include "TVector3.h"
#include "Nucleon.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"
#include "G4FermiPhaseSpaceDecay.hh"

using namespace std;

// Creating shortcut for an integer pair
typedef  pair<G4int, G4int> iPair;

// Structure to represent a graph
struct Graph
{
	// Vert and edges
    G4int V, E;
    vector< pair<G4double, iPair> > edges;

    // Adjacency lists
    vector< list<G4int> > adj;

    // Constructors
    Graph(G4int V, G4int E);
    Graph();

    // Destructor
    ~Graph();

    // Copy/Move
    Graph(const Graph& other);
    Graph& operator=(const Graph& other);
    Graph(Graph&& other) noexcept;
    Graph& operator=(Graph&& other) noexcept;

    // Utility function to add an edge
    void addEdge(G4int u, G4int v, G4double w);
    void DFSUtil(G4int v, G4bool visited[], vector < G4int > *clust_inside);
    // method to add a connection
    void addConn(G4int v, G4int w);
	// Method to print connected components in an
	// undirected graph
    vector< vector <G4int> > connectedComponents();

    // Function to find MST using Kruskal's
    // MST algorithm
    vector< vector <G4int> > AdvancedKruskalMST(G4double CD_in);
};

// To represent Disjoint Sets
struct DisjointSets
{
    G4int* parent, * rnk;
    G4int n;

    // Constructor
    DisjointSets(G4int n);

    // Destructor
    ~DisjointSets();

    // Find the parent of a node 'u'
    // Path Compression
    G4int find(G4int u);

    // Union by rank
    void merge(G4int x, G4int y);
};

class GMSTCluster{

	public:
	GMSTCluster(G4int Z_in, G4int A_in);
	~GMSTCluster();

	public:
	inline G4int GetZ() {return Z;};
	inline G4int GetA() {return A;};
	inline void SetZ(G4int Z_in) {Z = Z_in;}
	inline void SetA(G4int A_in) {A = A_in;}
	inline void PushBackCoordinateVector(TVector3 vec_in) {coord.push_back(vec_in);}
        inline std::vector<TVector3> GetCoordinates() {return coord;}

	private:
	G4int Z;
	G4int A;
	std::vector<TVector3> coord;
};

typedef std::vector<GMSTCluster> GMSTClusterVector;

class GMSTClustering{

	public:
	GMSTClustering();
	GMSTClustering(G4double CD_in, G4double Aa_in, G4double Ab_in);
	~GMSTClustering();

	public:
	inline void SetCD(G4double CD_in) {CritDist = CD_in; d0 = CD_in;}
	inline G4double GetCD(std::string side) {return (side == "A" ?  CritDistA : CritDist);}
    void SetCDExEn(G4double Ex, G4int A);
    void SetAlpha(G4double alpha_in){alphaPow = alpha_in;}
	Graph ClusterToGraph(aamcc::NucleonVector* nucleons, double V0pass);

	// ==== Phase-space MST controls ====
	inline void SetPhaseSpaceMode(bool on) { usePhaseSpace = on; }
	inline void SetV0(double V0_in)        { V0 = V0_in; }
	inline bool GetPhaseSpaceMode() const  { return usePhaseSpace; }
	inline double GetV0() const            { return V0; }

        // ==== Режим метрики ====
    // Условие по импульсу заменено условием по скорости: V = P/E (E — полная энергия).
    // Критерий DCM-подобный (per-particle), а не разность |V_i' - V_j'|:
    // 0 - чисто пространственный MST:
    //     w = Δr, ребро живо при w < dmin
    // 1 - метрика с per-particle cut в СЦМ пары (нуклон-нуклон), один проход при полном V0:
    //     w = max(Δr_eff, d_V_max · dmin / V0), ребро живо при w < dmin
    //     ⟺ Δr_eff < dmin И d_V_max < V0
    //     Δr_eff = √(Δx²+Δy²+(g·Δz)²), g — Лоренц-фактор пары из β_z=(pz_i+pz_j)/(E_i+E_j);
    //     d_V_max = max_k |V_k'|, V_k' = P_k'/E_k' — максимальная скорость нуклона пары
    //     в СЦМ пары (как в DCM BCOAV2 / Киреев T3)
    // 2 - та же метрика, но два иерархических прохода с постоянным V0 (см. GetClusters):
    //     1) порог V0, dmin = d0 (нуклон-нуклон, PairMetric);
    //     2) порог V0, dmin_eff = d0·(A_i+A_j)^(1/3) (кластер-кластер, PairMetric)
    //     PairMetric для кластеров проверяет все нуклоны объединения в его СЦМ
    //     (per-constituent). Вариант CrossCMMetric (проверка в СО другого кластера)
    //     реализован, но отвергнут: он строже и даёт ещё меньше тяжёлых фрагментов.
    inline void SetMetricMode(int m) { metricMode = m; }
    inline int  GetMetricMode() const { return metricMode; }

	std::vector<G4FragmentVector> GetClusters(aamcc::NucleonVector* nucleons, G4double ExA, G4double ExB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB);

	private:

    std::vector<G4FragmentVector> CalculateMomentum(std::vector<G4FragmentVector> noMomClusters, G4double ExEnA, G4double ExEnB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB, aamcc::NucleonVector rnucsA, aamcc::NucleonVector rnucsB, std::vector<int> rmapsA, std::vector<int> rmapsB);

    G4double CritDistA;
    G4double CritDist;
	G4double kappa = 1;
	G4double eps0 = 2.17*MeV;
	G4double alphaPow = -1.02;
	G4double d0 = -999;
	G4double aColRel = 5.315; //divided by aVol
	G4double aSurfRel  = 0.017; //divided by aVol
	G4double aVol     = 0.054;
    aamcc::NucleonVector* nucleonVector;

    G4FermiPhaseSpaceDecay phaseSpaceDecay;

    G4double SpecAa = 0;
    G4double SpecAb = 0;

    G4double a_opt = 2.243;
    G4double b_opt = 3.183 * MeV;
    G4double c_opt = 0.99;
    G4double d_opt = 0.29041;

    // ==== Phase-space MST controls (defaults) ====
    bool   usePhaseSpace = false;
    double V0 = 0.3;      // безразмерный порог |V'| в СЦМ объединения (per-particle, DCM-like)
    int metricMode = 1;   // 0 - spatial, 1 - pair-CM (1 pass), 2 - pair-CM (2 passes, constant V0)

    // Узел для метрики: нуклон (A=1) или кластер (A>1)
    struct Node {
        double x, y, z;
        double px, py, pz;
        double mass;  // полная масса покоя узла, МэВ
        int A;
    };

    // DCM-подобная метрика (per-particle cut в СЦМ объединения):
    // группа узлов a и группы b (нуклон или составные кластеры).
    // Ребро живо при Δr_eff < dmin И max_k |V_k'| < V0 в СЦМ объединения a+b.
    // Возвращает w = max(Δr_eff, d_V_max · dmin / V0).
    double PairMetric(const std::vector<Node>& a, const std::vector<Node>& b,
                      double dmin, double V0pass) const;

    // Альтернативная Cross-CM метрика (в mode 2 не используется):
    // нуклоны кластера a проверяются в СЦМ кластера b и наоборот
    // (кандидат проверяется в СО уже существующего кластера, как в DCM BCOAV2),
    // без "максимума по всем нуклонам объединения" (не подавляет тяжёлые остатки).
    double CrossCMMetric(const std::vector<Node>& a, const std::vector<Node>& b,
                         double dmin, double V0pass) const;

    // Второй проход: объединяет кластеры (списки индексов нуклонов),
    // используя СЦМ-координату/импульс кластеров и порог
    // dmin_eff = d0 · (A_i + A_j)^(1/3).
    std::vector<std::vector<G4int> > MergeClusters(
        const std::vector<std::vector<G4int> >& clusters,
        const aamcc::NucleonVector* nucleons,
        double d0, double V0pass) const;
};


#endif