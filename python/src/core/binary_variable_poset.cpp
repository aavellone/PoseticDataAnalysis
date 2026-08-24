
// MSVC: explicit standard includes (not pulled in transitively).
#include <vector>
#include <memory>
#include <utility>
#include <cstdint>
#include "my_exception.h"
#include "binary_variable_poset.h"
#include "linear_extension_generator_binary_variable.h"
#include "tensor.h"

#include <bit>
#include <format>
#include <string>
#include <algorithm>
#include <charconv>
#include <system_error>
#include <string_view>

// ***********************************************
// Factory Method
// ***********************************************

std::unique_ptr<BinaryVariablePOSet> BinaryVariablePOSet::Build(const std::vector<std::string>& variables) {
    if (variables.size() >= 64) {
        // Errore lanciato usando MyException e testo in inglese
        throw MyException("BinaryVariablePOSet supports a maximum of 63 variables.");
    }
    
    auto r = std::unique_ptr<BinaryVariablePOSet>(new BinaryVariablePOSet());
    r->variables_ = variables;
    r->numero_profili_ = 1ULL << r->variables_.size();
    
    r->element_names_.reserve(r->numero_profili_);
    
    for (std::uint64_t i = 0; i < r->numero_profili_; ++i) {
        r->element_names_.emplace_back(
                                    std::format("{:0{}b}", i, r->variables_.size())
                                    );
    }
    
    return r;
}

// ***********************************************
// Informazioni di Base
// ***********************************************

std::uint64_t BinaryVariablePOSet::NumberOfVariables() const noexcept {
    return variables_.size();
}

std::uint64_t BinaryVariablePOSet::size() const noexcept {
    return numero_profili_;
}

std::string BinaryVariablePOSet::to_string(char delimiter) const {
    return std::format("BinaryVariablePOSet(Variables: {}, Profiles: {})", variables_.size(), numero_profili_);
}

std::unique_ptr<LinearExtensionGenerator> BinaryVariablePOSet::CreateLinearExtensionGenerator()  {
    // LEGBinaryVariable richiede il numero di variabili (n_var).
    const std::uint64_t n_var = NumberOfVariables();
    return std::make_unique<LEGBinaryVariable>(n_var);
}

// ***********************************************
// Nomi e Identificatori
// ***********************************************

std::string_view BinaryVariablePOSet::GetElementName(std::uint64_t idx) const {
    if (idx >= numero_profili_) [[unlikely]] {
        throw MyException(std::format("Index {} out of bounds in GetElementName.", idx));
    }
    return element_names_[idx];
}

std::uint64_t BinaryVariablePOSet::GetElementId(std::string_view name) const  {
    if (name.empty()) [[unlikely]] {
        return 0ULL;
    }
    
    std::uint64_t result = 0;
    
    // std::from_chars è la funzione il parsing numerico.
    // L'ultimo parametro '2' indica la base binaria.
    auto [ptr, ec] = std::from_chars(name.data(), name.data() + name.size(), result, 2);
    
    // Controlliamo se c'è stato un errore (es. overflow o caratteri invalidi)
    // o se non tutta la stringa è stata consumata (es. "101abc" -> "abc" ignorato)
    if (ec != std::errc{} || ptr != name.data() + name.size()) [[unlikely]] {
        throw MyException(std::format("Invalid binary string '{}' in GetElementId.", name));
    }
    
    return result;
}


// ***********************************************
// Motore Relazionale On-The-Fly (Bitwise HPC)
// ***********************************************

bool BinaryVariablePOSet::IsLessOrEqual(std::uint64_t a, std::uint64_t b) const noexcept {
    return (a & b) == a; // a è un sottoinsieme di b
}

bool BinaryVariablePOSet::GreaterThan(std::uint64_t e1, std::uint64_t e2) const noexcept {
    return (e1 != e2) && ((e1 & e2) == e2); // e2 è un sottoinsieme proprio di e1
}


Tensor<std::uint8_t, 2> BinaryVariablePOSet::IncidenceMatrix() const {
    Tensor<std::uint8_t, 2> matrix({numero_profili_, numero_profili_}, 0);
    
    for (std::uint64_t sup = 0; sup < numero_profili_; ++sup) {
        matrix(sup, sup) = 1;
        
        for (std::uint64_t sub = (sup - 1) & sup; sub > 0; sub = (sub - 1) & sup) {
            matrix(sub, sup) = 1;
        }
        
        if (sup > 0) {
            matrix(0, sup) = 1;
        }
    }
    return matrix;
}


std::vector<std::pair<std::uint64_t, std::uint64_t>> BinaryVariablePOSet::Comparabilities() const {
    std::vector<std::pair<std::uint64_t, std::uint64_t>> comp;
    // Per un reticolo booleano le comparabilità totali (inclusi {0}) sono 3^V - 2^V
    // Omettiamo il reserve esatto per evitare calcoli con pow, affidandoci al vector.
    
    for (std::uint64_t sup = 1; sup < numero_profili_; ++sup) {
        // sub itererà in modo decrescente SOLO ed ESCLUSIVAMENTE sui sottoinsiemi di sup!
        for (std::uint64_t sub = (sup - 1) & sup; sub > 0; sub = (sub - 1) & sup) {
            comp.emplace_back(sub, sup);
        }
        // Lo zero è sottoinsieme di chiunque tranne che di se stesso (già escluso dal loop)
        comp.emplace_back(0, sup);
    }
    return comp;
}

// ***********************************************
// Utilità di Copia e Modifica
// ***********************************************

std::shared_ptr<POSet> BinaryVariablePOSet::Clone() const {
    auto c = std::unique_ptr<BinaryVariablePOSet>(new BinaryVariablePOSet());
    c->variables_ = this->variables_;
    c->numero_profili_ = this->numero_profili_;
    return c;
}


// ***********************************************
// ***********************************************
// ***********************************************

std::unique_ptr<POSet> BinaryVariablePOSet::Dual() const {
    auto dual = std::unique_ptr<BinaryVariablePOSet>(new BinaryVariablePOSet());
    // Il duale inverte la semantica dell'ordine: rovesciamo l'ordine delle variabili
    dual->variables_ = this->variables_;
    std::reverse(dual->variables_.begin(), dual->variables_.end());
    dual->numero_profili_ = this->numero_profili_;
    return dual;
}

// ***********************************************
// ***********************************************
// ***********************************************

[[nodiscard]] LatticeOfIdeals* BinaryVariablePOSet::GetLatticeOfIdeals() {
    throw MyException("GetLatticeOfIdeals() is not supported in BinaryVariablePOSet due to memory constraints. Use GetElementName() instead.");
}

// ***********************************************
// ***********************************************
// ***********************************************

[[nodiscard]] TreeOfIdeals* BinaryVariablePOSet::GetTreeOfIdeals() {
    throw MyException("GetTreeOfIdeals() is not supported in BinaryVariablePOSet due to memory constraints. Use GetElementName() instead.");
}

// ***********************************************
// ***********************************************
// ***********************************************

Tensor<std::uint8_t, 2> BinaryVariablePOSet::CoverMatrix() const {
    const std::uint64_t N = numero_profili_;
    
    Tensor<std::uint8_t, 2> result({N, N}, 0);
    
    // 3. Generazione bitwise diretta in O(V * 2^V) anziché O(N^2)
    for (std::uint64_t sub = 0; sub < N; ++sub) {
        
        // La maschera isola SOLO i bit a 0 di 'sub' (limitatamente alle V variabili attive).
        // Questi sono gli unici bit che, se accesi, generano un elemento che copre 'sub'.
        std::uint64_t mask = (~sub) & (N - 1);
        
        while (mask != 0) {
            // Estrae la posizione del bit meno significativo a 1
            int bit_pos = std::countr_zero(mask);
            
            // Calcola l'elemento che copre accendendo esattamente quel bit
            std::uint64_t sup = sub | (1ULL << bit_pos);
            
            // Imposta a 1 l'incidenza (riga = sub, colonna = sup)
            result(sub, sup) = 1;
            
            // Spegne il bit analizzato nella maschera per passare al prossimo
            mask &= mask - 1;
        }
    }
    
    return result;
}


// ***********************************************
// ***********************************************
// ***********************************************

bool BinaryVariablePOSet::IsTotalOrder() const {
    // Un reticolo booleano è una catena (ordine totale) solo se ha <= 1 variabile.
    return numero_profili_ <= 2;
}


// ***********************************************
// ***********************************************
// ***********************************************

bool BinaryVariablePOSet::IsExtensionOf(const POSet& p) const {
    // Usiamo dynamic_cast per verificare in modo sicuro se 'p'
    // è a tutti gli effetti un BinaryVariablePOSet.
    const auto* binary_p = dynamic_cast<const BinaryVariablePOSet*>(&p);
    
    // Se il cast ha successo (non è nullptr), allora 'p' è un BinaryVariablePOSet
    if (binary_p != nullptr) {
        // Due reticoli booleani sono un'estensione l'uno dell'altro
        // se e solo se sono isomorfi (cioè hanno lo stesso numero di variabili/profili).
        return this->numero_profili_ == binary_p->numero_profili_;
    }
    
    // Se 'p' è un POSet di tipo diverso, basandoci sulla tua logica originale,
    // assumiamo che non sia un'estensione diretta calcolabile in O(1).
    // Nota: Matematicamente bisognerebbe controllare che ogni arco di 'p'
    // esista in 'this', ma richiederebbe O(N^2) operazioni.
    throw MyException("IsExtensionOf is not supported for generic POSet types against BinaryVariablePOSet.");
}

// ***********************************************
// ***********************************************
// ***********************************************

BitSet BinaryVariablePOSet::DownSet(const std::vector<std::uint64_t>& els) const {
    // Allocazione singola e zero-initialized del BitSet dimensionato per il reticolo
    BitSet result(numero_profili_);
    
    // Iteriamo su tutti gli ID forniti in input
    for (std::uint64_t e : els) {
        // Guardia di sicurezza: controlliamo che l'elemento appartenga al reticolo
        if (e >= numero_profili_) {
            throw MyException(std::format("Element ID {} is out of bounds (max {}).", e, numero_profili_ - 1));
        }
        
        // Aggiungiamo l'elemento stesso al DownSet
        result.set(e);
        
        // Ciclo magico bitwise: itera *esclusivamente* sui sottoinsiemi di 'e'
        // in ordine decrescente, senza alcun "if" o branch.
        for (std::uint64_t sub = (e - 1) & e; sub > 0; sub = (sub - 1) & e) {
            // Il metodo set() del tuo BitSet gestisce automaticamente i duplicati
            // (se un sottoinsieme è già stato aggiunto da un altro elemento di 'els',
            //  l'operazione di OR bitwise lo lascerà invariato senza costi aggiuntivi).
            result.set(sub);
        }
        
        // L'insieme vuoto (0) è sempre nel DownSet di qualsiasi elemento (tranne 0 stesso)
        if (e > 0) {
            result.set(0);
        }
    }
    
    return result;
}

/**
 * @brief Checks if a given subset of elements forms a DownSet (Ideal).
 * * @details
 * **Mathematical Rule applied:** * In a Boolean lattice (isomorphic to a Power set lattice), a subset $S$ is a DownSet
 * if and only if for every element $e \in S$, all its direct lower covers (subsets
 * obtained by clearing exactly one active bit) are also present in $S$.
 * This avoids exponential inclusion-exclusion calculations, reducing complexity to $O(|S| \times V)$.
 * * **Bibliographic Reference:** * Davey, B.A., & Priestley, H.A. (2002). *Introduction to Lattices and Order* (2nd ed.).
 * Cambridge University Press. (See Chapter 1, Definition of Order Ideal / Down-set).
 * * @param els A vector containing the IDs of the elements in the subset.
 * @return true if the subset is a valid DownSet, false otherwise.
 * @throws MyException if any element ID is out of bounds.
 */
bool BinaryVariablePOSet::IsDownSet(const std::vector<std::uint64_t>& els) const {
    if (els.empty()) return true;
    
    BitSet s_bits(numero_profili_);
    for (std::uint64_t e : els) {
        if (e >= numero_profili_) {
            throw MyException(std::format("Element ID {} is out of bounds.", e));
        }
        s_bits.set(e);
    }
    
    for (std::uint64_t e : els) {
        std::uint64_t temp = e;
        while (temp != 0) {
            int bit_pos = std::countr_zero(temp);
            std::uint64_t child = e ^ (1ULL << bit_pos);
            if (!s_bits.test(child)) return false;
            temp &= temp - 1;
        }
    }
    return true;
}

/**
 * @brief Computes the UpSet (Filter) for a given subset of elements.
 * * @details
 * **Mathematical Rule applied:**
 * In a Boolean lattice, the UpSet of an element $e$ (principal filter) consists
 * of all its supersets. A superset is formed by taking $e$ and activating any
 * combination of the bits that are currently "off" in $e$. By isolating these
 * available bits into a `mask`, we can iterate exclusively over all its submasks
 * and use a bitwise OR to generate the exact supersets without any branch or
 * condition evaluation.
 * * **Bibliographic Reference:** * Davey, B.A., & Priestley, H.A. (2002). *Introduction to Lattices and Order* * (2nd ed.). Cambridge University Press. (See Chapter 1, Definition of Order Filter / Up-set).
 * * @param els A vector containing the IDs of the elements.
 * @return A BitSet where the active bits represent the elements in the UpSet.
 * @throws MyException if any element ID is out of bounds.
 */
BitSet BinaryVariablePOSet::UpSet(const std::vector<std::uint64_t>& els) const {
    BitSet result(numero_profili_);
    
    for (std::uint64_t e : els) {
        // Guardia di sicurezza
        if (e >= numero_profili_) {
            throw MyException(std::format("Element ID {} is out of bounds (max {}).", e, numero_profili_ - 1));
        }
        
        // 1. Isoliamo i bit che sono attualmente "spenti" in 'e'.
        // (numero_profili_ - 1) rappresenta la maschera di tutti '1' per le variabili attive (Universo).
        std::uint64_t mask = (~e) & (numero_profili_ - 1);
        
        // 2. Iteriamo esclusivamente su tutte le combinazioni possibili di questi bit disponibili.
        std::uint64_t sub = mask;
        do {
            // Generiamo il sovrainsieme combinando l'elemento base 'e' con la sottomaschera 'sub'
            result.set(e | sub);
            
            // Passiamo alla sottomaschera precedente.
            // Questa è un'istruzione bitwise nota in HPC per iterare le sottomaschere in ordine decrescente.
            sub = (sub - 1) & mask;
            
            // Quando sub diventa 0, l'istruzione (0 - 1) genera tutti 1 (underflow).
            // Il bitwise AND con 'mask' farà ritornare il valore a 'mask', interrompendo il ciclo.
        } while (sub != mask);
    }
    
    return result;
}

/**
 * @brief Retrieves the precomputed datastore of all UpSets (Not Supported).
 * * @details
 * **Mathematical Rule applied:**
 * In a Boolean lattice $\mathcal{P}(V)$, the total number of comparable pairs
 * (the sum of the sizes of all principal filters/UpSets) is exactly $3^V$.
 * Materializing this structure in memory requires exponential space $O(3^V)$,
 * which violates the $O(1)$ memory footprint of the on-the-fly bitwise architecture.
 * Therefore, this operation is explicitly disabled for this class.
 * * **Bibliographic Reference:**
 * Stanley, R. P. (2011). *Enumerative Combinatorics, Volume 1* (2nd ed.).
 * Cambridge University Press. (See Proposition 1.9.1 on Boolean algebras).
 * * @return A reference to the datastore (Will always throw).
 * @throws MyException Always, as full materialization is prohibited.
 * @note Ensure that `noexcept` is removed from the base `POSet` interface
 * so this exception can propagate safely.
 */
const POSet::DATASTORE& BinaryVariablePOSet::UpSets() const  {
    throw MyException("Precomputing all UpSets is not supported in BinaryVariablePOSet due to O(3^V) memory complexity. Use on-the-fly bitwise methods instead.");
}

/**
 * @brief Checks if a given subset of elements forms an UpSet (Filter).
 *
 * @details
 * **Mathematical Rule applied:**
 * In a Boolean lattice (isomorphic to a Power set lattice), a subset $S$ is an UpSet
 * if and only if for every element $e \in S$, all its direct upper covers (supersets
 * obtained by activating exactly one available bit) are also present in $S$.
 * This topological check avoids exponential inclusion-exclusion algorithms,
 * bringing the time complexity down to $O(|S| \times V)$.
 *
 * **Bibliographic Reference:**
 * Davey, B.A., & Priestley, H.A. (2002). *Introduction to Lattices and Order* * (2nd ed.). Cambridge University Press. (See Chapter 1, Definition of Order Filter / Up-set).
 *
 * @param els A vector containing the IDs of the elements in the subset.
 * @return true if the subset is a valid UpSet, false otherwise.
 * @throws MyException if any element ID is out of bounds.
 */
bool BinaryVariablePOSet::IsUpSet(const std::vector<std::uint64_t>& els) const {
    if (els.empty()) {
        return true; // Per definizione, l'insieme vuoto è un UpSet
    }
    
    // 1. Carichiamo tutti gli elementi in un BitSet per un lookup in O(1)
    BitSet s_bits(numero_profili_);
    for (std::uint64_t e : els) {
        if (e >= numero_profili_) {
            throw MyException(std::format("Element ID {} is out of bounds (max {}).", e, numero_profili_ - 1));
        }
        s_bits.set(e);
    }
    
    // Maschera dell'Universo: tutti i bit corrispondenti alle variabili sono a 1
    const std::uint64_t universe_mask = numero_profili_ - 1;
    
    // 2. Verifica topologica dell'UpSet (Zero branch-misprediction)
    for (std::uint64_t e : els) {
        // Isoliamo ESCLUSIVAMENTE i bit che sono attualmente SPENTI nell'elemento 'e'
        std::uint64_t available_bits = (~e) & universe_mask;
        
        // Iteriamo solo sui bit disponibili
        while (available_bits != 0) {
            // Troviamo la posizione del bit meno significativo spento
            int bit_pos = std::countr_zero(available_bits);
            
            // Generiamo il sovrainsieme diretto (parent) ACCENDENDO quel bit
            std::uint64_t parent = e | (1ULL << bit_pos);
            
            // Se un sovrainsieme diretto non fa parte dell'insieme,
            // la proprietà di chiusura verso l'alto è violata.
            if (!s_bits.test(parent)) {
                return false;
            }
            
            // Rimuoviamo il bit visitato per passare al prossimo
            available_bits &= available_bits - 1;
        }
    }
    
    // Se non ci sono state violazioni, l'insieme è un UpSet valido.
    return true;
}


/**
 * @brief Generates the complete set of order relations (comparabilities) as string pairs.
 *
 * @details
 * **Mathematical Rule applied:**
 * The total number of comparability pairs $(a, b)$ where $a \le b$ in a Boolean lattice
 * $\mathcal{P}(V)$ is exactly $3^V$. The old $O(N^2)$ algorithm checked all $4^V$ possible
 * pairs. This HPC implementation iterates through each element and exclusively visits its
 * true subsets using a bitwise descending mask, achieving optimal $O(3^V)$ time complexity.
 * Due to the exponential memory required to materialize string pairs, a safety guard
 * is enforced for $V > 15$.
 *
 * **Bibliographic Reference:**
 * Stanley, R. P. (2011). *Enumerative Combinatorics, Volume 1* (2nd ed.).
 * Cambridge University Press. (See Section 1.9 on Boolean algebras and $3^V$ intervals).
 *
 * @return A vector of string pairs representing all $(a, b)$ relations.
 * @throws MyException if variables > 15, to prevent Out-Of-Memory crashes.
 */
std::vector<std::pair<std::string, std::string>> BinaryVariablePOSet::OrderRelation() const {
    std::uint64_t V = variables_.size();
    
    
    std::uint64_t exact_capacity = 1;
    for (std::uint64_t i = 0; i < V; ++i) {
        exact_capacity *= 3;
    }
    
    std::vector<std::pair<std::string, std::string>> result;
    result.reserve(exact_capacity);
    
    for (std::uint64_t sup = 0; sup < numero_profili_; ++sup) {
        std::string_view sup_name = GetElementName(sup);
        
        // Aggiungiamo la relazione riflessiva (a <= a)
        result.emplace_back(std::string(sup_name), std::string(sup_name));
        
        // Iteriamo magicamente in discesa SOLO sui veri sottoinsiemi di 'sup'
        for (std::uint64_t sub = (sup - 1) & sup; sub > 0; sub = (sub - 1) & sup) {
            result.emplace_back(std::string(GetElementName(sub)), std::string(sup_name));
        }
        
        // Gestione speciale dello zero (sottoinsieme di tutti tranne che di se stesso)
        if (sup > 0) {
            result.emplace_back(std::string(GetElementName(0)), std::string(sup_name));
        }
    }
    
    return result;
}

/**
 * @brief Computes the Cover Relation (Hasse diagram edges) of the Boolean lattice.
 *
 * @details
 * **Mathematical Rule applied:**
 * In a Boolean lattice $\mathcal{P}(V)$, an element $A$ is covered by $B$ ($A \prec B$)
 * if and only if $A \subset B$ and $|B \setminus A| = 1$. The total number of cover
 * relations is exactly $V \cdot 2^{V-1}$. The optimized algorithm avoids the $O(N^2)$
 * pair combinations by exclusively iterating over the unset bits of each element
 * and generating its direct upper covers.
 * Note: The mathematical definition of a cover relation is irreflexive, hence
 * reflexive pairs (a, a) are strictly omitted.
 *
 * **Bibliographic Reference:**
 * Davey, B.A., & Priestley, H.A. (2002). *Introduction to Lattices and Order* * (2nd ed.). Cambridge University Press. (See Chapter 1, Definition of Covering).
 *
 * @return A vector of pairs (subset, superset) representing the cover relations.
 * @throws MyException if the memory required exceeds safe limits (V > 24).
 */
std::vector<std::pair<std::uint64_t, std::uint64_t>> BinaryVariablePOSet::CoverRelation() const {
    std::uint64_t V = variables_.size();
    
    // 1. Guardia di sicurezza per la memoria
    // Ogni std::pair di uint64_t occupa 16 byte.
    // Per V=24, avremmo 24 * 2^23 = 201.326.592 archi, che occupano circa 3.2 GB di RAM.
    if (V > 24) {
        throw MyException(std::format("Cannot materialize CoverRelation for {} variables: exceeded safe memory limits.", V));
    }
    
    // 2. Calcolo esatto della capacità: V * 2^(V-1)
    std::uint64_t exact_capacity = V * (numero_profili_ >> 1);
    
    std::vector<std::pair<std::uint64_t, std::uint64_t>> result;
    result.reserve(exact_capacity); // Pre-allocazione per zero memory-fragmentation
    
    // Maschera dell'Universo
    const std::uint64_t universe_mask = numero_profili_ - 1;
    
    // 3. Generazione diretta degli archi del Diagramma di Hasse
    for (std::uint64_t sub = 0; sub < numero_profili_; ++sub) {
        // Isoliamo ESCLUSIVAMENTE i bit SPENTI in 'sub'
        std::uint64_t available_bits = (~sub) & universe_mask;
        
        while (available_bits != 0) {
            // Troviamo il primo bit spento (hardware nativo O(1))
            int bit_pos = std::countr_zero(available_bits);
            
            // Generiamo l'elemento che copre 'sub' ACCENDENDO quel bit
            std::uint64_t sup = sub | (1ULL << bit_pos);
            
            // Inseriamo la relazione (sottoinsieme, sovrainsieme)
            result.emplace_back(sub, sup);
            
            // Rimuoviamo il bit visitato per iterare al successivo
            available_bits &= available_bits - 1;
        }
    }
    
    return result;
}

/**
 * @brief Computes the Comparability Set of a single element.
 *
 * @details
 * **Mathematical Rule applied:**
 * In a partially ordered set, the comparability set of an element $e$ is defined
 * as $\{x \in P \mid x \le e \lor e \le x\}$. In our Boolean lattice $\mathcal{P}(V)$,
 * this is the direct union of its principal ideal (DownSet, its subsets) and its
 * principal filter (UpSet, its supersets).
 * This HPC implementation bypasses container allocations entirely, generating both
 * sets on-the-fly via bitwise descending masks and natively storing the union in
 * the resulting BitSet.
 *
 * **Bibliographic Reference:**
 * Davey, B.A., & Priestley, H.A. (2002). *Introduction to Lattices and Order* * (2nd ed.). Cambridge University Press. (See Chapter 1, Definition of Comparability).
 *
 * @param e The ID of the element to analyze.
 * @return A BitSet where active bits represent elements comparable to $e$.
 * @throws MyException if the element ID is out of bounds.
 */
BitSet BinaryVariablePOSet::ComparabilitySetOf(std::uint64_t e) const {
    // Guardia di sicurezza
    if (e >= numero_profili_) {
        throw MyException(std::format("Element ID {} is out of bounds (max {}).", e, numero_profili_ - 1));
    }
    
    BitSet result(numero_profili_);
    
    // 1. Generazione del DownSet (tutti i sottoinsiemi di 'e')
    // Iteriamo magicamente in discesa usando (sub - 1) & e
    std::uint64_t sub = e;
    do {
        result.set(sub);
        sub = (sub - 1) & e;
        // Quando sub diventa 0, (0 - 1) genera underflow (tutti 1).
        // Il bitwise AND con 'e' riporta il valore a 'e', terminando il ciclo.
    } while (sub != e);
    
    // 2. Generazione dell'UpSet (tutti i sovrainsiemi di 'e')
    // Isoliamo i bit che sono attualmente "spenti" in 'e'
    std::uint64_t mask = (~e) & (numero_profili_ - 1);
    std::uint64_t sup = mask;
    do {
        result.set(e | sup);
        sup = (sup - 1) & mask;
    } while (sup != mask); // Stessa logica di underflow per la condizione di uscita
    
    return result;
}

/**
 * @brief Computes the Incomparability Set of a single element using global complement.
 *
 * @details
 * **Mathematical Rule applied:**
 * In any partially ordered set, the Incomparability Set of an element $e$ is exactly
 * the complement of its Comparability Set: $P \setminus (\{x \in P \mid x \le e \lor e \le x\})$.
 * Since the BitSet class provides an HPC-optimized flip() method that performs
 * an in-place bitwise NOT (in O(N/64)) and safely masks out-of-bounds bits,
 * we can reuse ComparabilitySetOf to achieve massive performance gains over
 * the traditional O(N) iterative approach.
 *
 * **Bibliographic Reference:**
 * Davey, B.A., & Priestley, H.A. (2002). *Introduction to Lattices and Order* (2nd ed.).
 * Cambridge University Press. (See Chapter 1, relations and complements).
 *
 * @param e The ID of the element to analyze.
 * @return A BitSet where active bits represent elements incomparable to e.
 * @throws MyException if the element ID is out of bounds.
 */
BitSet BinaryVariablePOSet::IncomparabilitySetOf(std::uint64_t e) const {
    if (e >= numero_profili_) {
        throw MyException(std::format("Element ID {} is out of bounds (max {}).", e, numero_profili_ - 1));
    }
    
    // 1. Calcola l'insieme di comparabilità (DownSet U UpSet) in tempo O(C)
    // C è il numero di elementi comparabili, che è infinitesimo rispetto a N.
    BitSet result = ComparabilitySetOf(e);
    
    // 2. Inverte tutti i bit dell'array flat in tempo vettoriale O(N/64).
    // Il metodo flip() del tuo BitSet ri-azzera automaticamente i bit di padding,
    // mantenendo il risultato matematicamente impeccabile.
    result.flip();
    
    return result;
}

/**
 * @brief Retrieves the set of maximal elements in the lattice.
 *
 * @details
 * **Mathematical Rule applied:**
 * A Boolean lattice $\mathcal{P}(V)$ is a bounded lattice. It possesses a unique
 * greatest element (the Top element, $\top$), which corresponds to the universe set
 * containing all variables. Its integer representation is strictly $2^V - 1$.
 * Therefore, the set of maximals contains exactly one element, retrievable in $O(1)$.
 *
 * **Bibliographic Reference:**
 * Davey, B.A., & Priestley, H.A. (2002). *Introduction to Lattices and Order* (2nd ed.).
 * Cambridge University Press. (See Chapter 2, Bounded Lattices).
 *
 * @return A BitSet containing only the maximal element (Top).
 */
BitSet BinaryVariablePOSet::Maximals() const {
    BitSet result(numero_profili_);
    // Imposta a 1 solo l'ultimo ID del reticolo (L'Universo)
    result.set(numero_profili_ - 1);
    return result;
}

/**
 * @brief Retrieves the set of minimal elements in the lattice.
 *
 * @details
 * **Mathematical Rule applied:**
 * Similar to the Top element, a Boolean lattice possesses a unique least element
 * (the Bottom element, $\bot$), which corresponds to the empty set.
 * Its integer representation is strictly $0$. The set of minimals is simply $\{0\}$.
 *
 * **Bibliographic Reference:**
 * Davey, B.A., & Priestley, H.A. (2002). *Introduction to Lattices and Order* (2nd ed.).
 * Cambridge University Press. (See Chapter 2, Bounded Lattices).
 *
 * @return A BitSet containing only the minimal element (Bottom).
 */
BitSet BinaryVariablePOSet::Minimals() const {
    BitSet result(numero_profili_);
    // Imposta a 1 solo lo zero (L'Insieme Vuoto)
    result.set(0);
    return result;
}

/**
 * @brief Checks if a given element is a maximal element.
 *
 * @param e The ID of the element to check.
 * @return true if the element is the maximal element (Top), false otherwise.
 * @throws MyException if the element ID is out of bounds.
 */
bool BinaryVariablePOSet::IsMaximal(std::uint64_t e) const {
    if (e >= numero_profili_) {
        throw MyException(std::format("Element ID {} is out of bounds (max {}).", e, numero_profili_ - 1));
    }
    return e == (numero_profili_ - 1);
}

/**
 * @brief Checks if a given element is a minimal element.
 *
 * @param e The ID of the element to check.
 * @return true if the element is the minimal element (Bottom), false otherwise.
 * @throws MyException if the element ID is out of bounds.
 */
bool BinaryVariablePOSet::IsMinimal(std::uint64_t e) const {
    if (e >= numero_profili_) {
        throw MyException(std::format("Element ID {} is out of bounds (max {}).", e, numero_profili_ - 1));
    }
    return e == 0;
}

/**
 * @brief Computes all pairs of incomparable elements in the Boolean lattice.
 *
 * @details
 * **Mathematical Rule applied:**
 * Two elements $a$ and $b$ are incomparable ($a \parallel b$) if $a \not\le b$ and $b \not\le a$.
 * In a Boolean lattice $\mathcal{P}(V)$, generating incomparable pairs is highly dense.
 * The total number of undirected incomparable pairs is given by the exact formula:
 * $\frac{4^V - 2 \cdot 3^V + 2^V}{2}$.
 * By taking advantage of numerical ordering in the loop ($v_1 < v_2$), we guarantee
 * $v_2 \not\subseteq v_1$ a priori. Thus, we only need a single bitwise check
 * (`(v1 & ~v2) != 0`) to verify $v_1 \not\subseteq v_2$, cutting the evaluations in half.
 *
 * **Bibliographic Reference:**
 * Stanley, R. P. (2011). *Enumerative Combinatorics, Volume 1* (2nd ed.).
 * Cambridge University Press. (See combinatorial counting on Boolean algebras).
 *
 * @return A vector of pairs representing all incomparable relationships (v1, v2) where v1 < v2.
 * @throws MyException if the variables exceed V=14 due to extreme memory constraints.
 */
std::vector<std::pair<std::uint64_t, std::uint64_t>> BinaryVariablePOSet::Incomparabilities() const {
    std::uint64_t V = variables_.size();
    
    // 1. Guardia di sicurezza per la memoria (Evita l'OOM Killer)
    // Per V=14, il vector allocherà circa 2 GB.
    // Per V=15, balzeremmo a quasi 8 GB.
    if (V > 14) {
        throw MyException(std::format("Cannot materialize Incomparabilities for {} variables: exceeded maximum memory limit (max V=14).", V));
    }
    
    // 2. Pre-calcolo esatto della capacità per azzerare le riallocazioni
    std::uint64_t pow3 = 1;
    for (std::uint64_t i = 0; i < V; ++i) {
        pow3 *= 3;
    }
    
    std::uint64_t pow4 = 1ULL << (2 * V); // 4^V
    std::uint64_t pow2 = 1ULL << V;       // 2^V
    
    // Formula per le coppie non orientate incomparabili
    std::uint64_t exact_capacity = (pow4 - 2 * pow3 + pow2) / 2;
    
    std::vector<std::pair<std::uint64_t, std::uint64_t>> result;
    result.reserve(exact_capacity);
    
    // 3. Generazione Bruteforce ma ottimizzata al singolo ciclo macchina
    for (std::uint64_t v1 = 0; v1 < numero_profili_; ++v1) {
        // v1 < v2 garantisce implicitamente che v2 non possa essere sottoinsieme di v1.
        for (std::uint64_t v2 = v1 + 1; v2 < numero_profili_; ++v2) {
            
            // Singolo check branchless: v1 possiede almeno un bit che v2 non ha?
            if ((v1 & ~v2) != 0) {
                result.emplace_back(v1, v2);
            }
            
        }
    }
    
    return result;
}

/**
 * @brief Generates the first canonical Linear Extension of the Boolean lattice.
 *
 * @details
 * **Mathematical Rule applied:**
 * A Linear Extension is a total ordering $L = (x_1, x_2, \dots, x_n)$ such that
 * $x_i \le_{\text{poset}} x_j \implies i \le j$. In a Boolean lattice where elements
 * are identified by their bitmask, the natural numerical order of integers
 * is a valid linear extension. This is because $A \subseteq B \implies \text{int}(A) \le \text{int}(B)$.
 *
 * **Bibliographic Reference:**
 * Stanley, R. P. (2011). *Enumerative Combinatorics, Volume 1* (2nd ed.).
 * Cambridge University Press. (See Chapter 3, Linear Extensions of Posets).
 *
 * @param[out] le The LinearExtension object to be populated.
 */
void BinaryVariablePOSet::FirstLE(LinearExtension& le) const {
    const std::uint64_t n = le.size();
    
    for (std::uint64_t k = 0; k < n; ++k) {
        le.Set(k, k);
    }
}

/**
 * @brief Checks if two elements are comparable in the Boolean lattice.
 *
 * @details
 * **Mathematical Rule applied:**
 * In the Boolean lattice $\mathcal{P}(V)$, two elements $a, b$ are comparable
 * ($a \sim b$) if and only if $a \le b$ or $b \le a$.
 * This is equivalent to $a \subseteq b \lor b \subseteq a$.
 * Using bitwise logic, this is efficiently checked as:
 * `((a & ~b) == 0) || ((b & ~a) == 0)`.
 *
 * **Bibliographic Reference:**
 * Davey, B.A., & Priestley, H.A. (2002). *Introduction to Lattices and Order* (2nd ed.).
 * Cambridge University Press. (See Chapter 1, Basic Definitions).
 *
 * @param a ID of the first element.
 * @param b ID of the second element.
 * @return true if a and b are comparable, false otherwise.
 */
[[nodiscard]] bool BinaryVariablePOSet::IsComparable(std::uint64_t a, std::uint64_t b) const noexcept {
    // 1. Verifica se a è sottoinsieme di b: (a & ~b) == 0
    // 2. Verifica se b è sottoinsieme di a: (b & ~a) == 0
    // L'operatore || in C++ è a corto circuito: se la prima condizione è vera,
    // la seconda non viene nemmeno valutata dalla CPU.
    
    return ((a & ~b) == 0) || ((b & ~a) == 0);
}

/**
 * @brief Computes the probabilistic dominance matrix using the LPOM (Relative Ideals) model.
 * * @details This implementation uses the original Local Partial Order Model formula
 * (known in literature as Equation 13), developed by Brüggemann, Lerche, and Sørensen.
 * For pairs of incomparable elements, the mutual probability is estimated by evaluating
 * the cardinalities of the relative ideals, meaning the set difference between the
 * Up-sets and Down-sets of the two elements: |U(x) \ U(y)| and |O(x) \ O(y)|.
 * * In the context of a BinaryVariablePOSet (Boolean Lattice), set operations are
 * hardware-accelerated using bitwise arithmetic and SIMD instructions (e.g., popcount).
 * * @note From a computational standpoint, this variant requires calculating the logical
 * intersection for every single pair (v1, v2) within the innermost loop, resulting in
 * an O(N^2) complexity for bitwise intersections.
 * * @see NERI Technical Report No. 479 (2004) - Equation 13.
 * @see Brüggemann, R., Sørensen, P. B., Lerche, D. B., Carlsen, L. (2004).
 * "Estimation of Averaged Ranks by a Local Partial Order Model".
 * * @return Tensor<double, 2> A dense square matrix where element (i, j) represents
 * the probabilistic estimate P(i > j) that element 'i' dominates element 'j'.
 */
[[nodiscard]] Tensor<double, 2> BinaryVariablePOSet::BLSDominanceRelative() const {
    const std::uint64_t num_vars = variables_.size();

    // 2. Allocazione locale della matrice (ritorno per valore)
    Tensor<double, 2> matrice({numero_profili_, numero_profili_}, kUninitialized);
    
    // 3. Imposta la diagonale (ogni elemento domina se stesso al 100%)
    for (std::uint64_t k = 0; k < numero_profili_; ++k) {
        matrice(k, k) = 1.0;
    }
    
    // 4. Calcolo della dominanza incrociata
    for (std::uint64_t v1 = 0; v1 < numero_profili_; ++v1) {
        // Pre-calcoliamo i bit attivi di v1 (operazione hardware nativa)
        const std::uint64_t pop_v1 = static_cast<std::uint64_t>(std::popcount(v1));
        
        for (std::uint64_t v2 = v1 + 1; v2 < numero_profili_; ++v2) {
            const std::uint64_t v_or = (v1 | v2);
            
            if (v_or == v2) {
                // v1 è un sottoinsieme di v2 -> dominanza diretta
                matrice(v1, v2) = 1.0;
                matrice(v2, v1) = 0.0;
            } else {
                // Elementi Incomparabili: Formula BLS per Reticolo Booleano (Powerset)
                const std::uint64_t v_and = (v1 & v2);
                
                const std::uint64_t pop_v2  = static_cast<std::uint64_t>(std::popcount(v2));
                const std::uint64_t pop_or  = static_cast<std::uint64_t>(std::popcount(v_or));
                const std::uint64_t pop_and = static_cast<std::uint64_t>(std::popcount(v_and));
                // Calcolo upij = |U(v1) \ U(v2) \ {v1, v2}|
                // U(x) ha dimensione 2^(V - popcount(x)).
                // L'intersezione U(v1) e U(v2) è l'Up-set della loro unione logica (v_or).
                // Togliamo 1 per escludere v1 stesso.
                const double upij = static_cast<double>((1ULL << (num_vars - pop_v1)) - (1ULL << (num_vars - pop_or)) - 1);
                const double upji = static_cast<double>((1ULL << (num_vars - pop_v2)) - (1ULL << (num_vars - pop_or)) - 1);
                
                // Calcolo downij = |O(v1) \ O(v2) \ {v1, v2}|
                // O(x) ha dimensione 2^popcount(x).
                // L'intersezione O(v1) e O(v2) è il Down-set della loro intersezione logica (v_and).
                // Togliamo 1 per escludere v1 stesso.
                const double downij = static_cast<double>((1ULL << pop_v1) - (1ULL << pop_and) - 1);
                const double downji = static_cast<double>((1ULL << pop_v2) - (1ULL << pop_and) - 1);
                
                // Rapporti (il +1.0 della formula originale)
                const double aij = (upij + 1.0) / (downij + 1.0);
                const double aji = (upji + 1.0) / (downji + 1.0);
                
                const double sum_a = aij + aji;
                
                // Assegnazione simmetrica incrociata sulla matrice
                matrice(v2, v1) = aji / sum_a;
                matrice(v1, v2) = aij / sum_a;
            }
        }
    }
    
    return matrice;
}

/**
 * @brief Computes the probabilistic dominance matrix using the LPOM (Absolute Ideals) model.
 * * @details Implements Equation 13' (the Q'-based variant) from the theory by Brüggemann et al.
 * Unlike the relative variant, this function estimates the degrees of freedom of an element
 * based exclusively on the cardinality of its absolute ideals: |U(x)| and |O(x)|.
 * * @note HPC Optimization: This metric offers vastly superior performance because the model
 * parameters' dependency is decoupled. This allows the extraction of the calculations for
 * element 'v1' from the innermost nested loop, reducing the heavy floating-point operations
 * from O(N^2) to O(N) and completely eliminating the need to evaluate cross-set intersections.
 * * @see NERI Technical Report No. 479 (2004) - Equation 13'.
 * * @return Tensor<double, 2> A dense square matrix where element (i, j) represents
 * the probabilistic estimate P(i > j) calculated using absolute ideals.
 */
[[nodiscard]] Tensor<double, 2> BinaryVariablePOSet::BLSDominanceAbsolute() const {
    const std::uint64_t num_vars = variables_.size();
    
    Tensor<double, 2> matrice({numero_profili_, numero_profili_}, kUninitialized);
    
    for (std::uint64_t k = 0; k < numero_profili_; ++k) {
        matrice(k, k) = 1.0;
    }
    
    for (std::uint64_t v1 = 0; v1 < numero_profili_; ++v1) {
        const std::uint64_t pop_v1 = static_cast<std::uint64_t>(std::popcount(v1));
        
        // Calcolo ideali assoluti per v1: |U(v1)| e |O(v1)| (meno se stesso)
        const double up_size_v1 = static_cast<double>((1ULL << (num_vars - pop_v1)) - 1);
        const double down_size_v1 = static_cast<double>((1ULL << pop_v1) - 1);
        
        // Q' per v1
        const double q_v1 = (up_size_v1 + 1.0) / (down_size_v1 + 1.0);
        
        for (std::uint64_t v2 = v1 + 1; v2 < numero_profili_; ++v2) {
            const std::uint64_t v_or = (v1 | v2);
            
            if (v_or == v2) {
                matrice(v1, v2) = 1.0;
                matrice(v2, v1) = 0.0;
            } else {
                const std::uint64_t pop_v2 = static_cast<std::uint64_t>(std::popcount(v2));
                
                // Calcolo ideali assoluti per v2
                const double up_size_v2 = static_cast<double>((1ULL << (num_vars - pop_v2)) - 1);
                const double down_size_v2 = static_cast<double>((1ULL << pop_v2) - 1);
                
                // Q' per v2
                const double q_v2 = (up_size_v2 + 1.0) / (down_size_v2 + 1.0);
                
                const double sum_q = q_v1 + q_v2;
                
                // Assegnazione matrice Eq 13' (P(v1 > v2) = q_v2 / (q_v1 + q_v2))
                // riga dominata dalla colonna
                matrice(v2, v1) = q_v2 / sum_q;
                matrice(v1, v2) = q_v1 / sum_q;
            }
        }
    }
    
    return matrice;
}
