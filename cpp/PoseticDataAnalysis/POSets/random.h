#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <utility>
#include <cstddef>

// =============================================================================
// random.h — generatore pseudo-casuale del progetto.
//
// Motore: SplitMix64 con aritmetica di bit esplicita. Sostituisce il precedente
// std::mt19937_64 + std::uniform_int/real_distribution + std::shuffle: quelle
// distribuzioni NON sono portabili tra librerie standard (libstdc++ e libc++
// producono sequenze diverse a parita' di engine e seed). Poiche' il core viene
// distribuito come package R e Python compilati su piu' piattaforme, la
// riproducibilita' cross-platform e' essenziale: SplitMix64 con bit-math e'
// bit-identico ovunque, minuscolo e senza dipendenze dal <random> sul percorso
// caldo.
//
// L'INTERFACCIA PUBBLICA e' invariata rispetto alla versione mt19937 (stessi
// nomi di metodo usati dal core: InitBubleyDyer / RndUpdate / RndPosition /
// RndNext / RndNextInt / RndShuffle / RndSample / Restart / Seed / GLOBAL /
// STARTUP_SEED). Cambia solo lo stream prodotto, che ora e' portabile.
// =============================================================================

#include <cstdint>
#include <vector>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <limits>

#if defined(_MSC_VER)
#include <intrin.h>  // __umulh: MSVC non ha __uint128_t
#endif

#include "my_exception.h"

class Random {
private:
    uint64_t _seed;
    uint64_t _state;              // stato SplitMix64
    uint64_t _positionRange = 0;  // RndPosition estrae in [0, _positionRange]

    // Buffer di bit per RndUpdate: una singola estrazione (64 bit) fornisce 64
    // coin flip, riducendo di ~64x le chiamate al generatore.
    uint64_t _bitBuffer = 0;
    uint32_t _bitsLeft = 0;

    static constexpr uint64_t kGamma = 0x9E3779B97F4A7C15ULL;  // incremento (golden ratio)

    /// Passo di avalanche finale di SplitMix64.
    static inline uint64_t Finalize(uint64_t z) noexcept {
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    /// Prossimo valore a 64 bit (un passo SplitMix64; avanza lo stream).
    inline uint64_t NextU64() noexcept {
        _state += kGamma;
        return Finalize(_state);
    }

    /// Intero uniforme in [0, n) via moltiplicazione-shift di Lemire (bias di
    /// modulo trascurabile). Ritorna 0 se n == 0.
    /// Calcola i 64 bit alti del prodotto a 64x64 bit: identico su tutte le
    /// piattaforme (GCC/Clang via __uint128_t, MSVC via __umulh).
    inline uint64_t Below(uint64_t n) noexcept {
        if (n == 0) return 0;
#if defined(_MSC_VER)
        return __umulh(NextU64(), n);
#else
        return static_cast<uint64_t>(
            (static_cast<__uint128_t>(NextU64()) * static_cast<__uint128_t>(n)) >> 64);
#endif
    }

public:
    inline static const uint64_t STARTUP_SEED = static_cast<uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());

    static Random GLOBAL;

    explicit Random(uint64_t s) : _seed(s), _state(s) {}
    Random() : _seed(STARTUP_SEED), _state(STARTUP_SEED) {}

    Random(const Random&)            = delete;
    Random& operator=(const Random&) = delete;

    Random(Random&&) noexcept            = default;
    Random& operator=(Random&&) noexcept = default;

    /// Configura l'intervallo di RndPosition su [0, size-2] (nessun consumo di
    /// stream): per Bubley-Dyer, la posizione dello swap adiacente.
    inline void InitBubleyDyer(uint64_t size) {
        _positionRange = (size > 1) ? (size - 2) : 0;
    }

    /// Coin flip (0/1); consuma 64 bit da una sola estrazione, LSB per primo.
    inline uint64_t RndUpdate() {
        if (_bitsLeft == 0) {
            _bitBuffer = NextU64();
            _bitsLeft = 64;
        }
        const uint64_t bit = _bitBuffer & 1ULL;
        _bitBuffer >>= 1;
        --_bitsLeft;
        return bit;
    }

    /// Posizione uniforme in [0, size-2] (impostata da InitBubleyDyer).
    inline uint64_t RndPosition() {
        return Below(_positionRange + 1);
    }

    [[nodiscard]] inline uint64_t Seed() const noexcept { return _seed; }

    /// Riporta lo stato al seed iniziale e svuota il buffer dei coin flip: la
    /// stessa sequenza viene riprodotta dopo il riavvio (stato riproducibile).
    inline void Restart() noexcept {
        _state = _seed;
        _bitBuffer = 0;
        _bitsLeft = 0;
    }

    /// Double uniforme in [0, 1) con risoluzione a 53 bit di mantissa.
    [[nodiscard]] inline double RndNext() {
        return static_cast<double>(NextU64() >> 11) * (1.0 / 9007199254740992.0);
    }

    /// Intero uniforme nell'intervallo CHIUSO [min, max] (precondizione min<=max).
    [[nodiscard]] inline uint64_t RndNextInt(uint64_t min, uint64_t max) {
        const uint64_t span = max - min;
        // Range pieno a 64 bit: span+1 andrebbe in overflow a 0, quindi si usa
        // direttamente un'estrazione a 64 bit (uniforme su tutto l'intervallo).
        if (span == std::numeric_limits<uint64_t>::max()) {
            return NextU64();
        }
        return min + Below(span + 1);
    }

    /// Fisher-Yates in-place su un vettore di uint64.
    inline void RndShuffle(std::vector<uint64_t>& v) {
        for (std::size_t i = v.size(); i > 1; --i) {
            const std::size_t j = static_cast<std::size_t>(Below(static_cast<uint64_t>(i)));
            std::swap(v[i - 1], v[j]);
        }
    }

    /// Campiona senza reinserimento k interi da [min, max].
    [[nodiscard]] inline std::vector<uint64_t> RndSample(uint64_t min, uint64_t max, uint64_t k = 1) {
        if (min > max) throw MyException("Random::RndSample: min > max");

        const uint64_t pop = max - min + 1;
        const uint64_t n   = std::min(k, pop);

        std::vector<uint64_t> pool(pop);
        std::iota(pool.begin(), pool.end(), min);

        for (uint64_t i = 0; i < n; ++i) {
            const uint64_t pick = i + Below(pop - i);  // [i, pop-1]
            std::swap(pool[i], pool[static_cast<std::size_t>(pick)]);
        }

        pool.resize(n);
        return pool;
    }
};

// ==========================================================
// DEFINIZIONE VARIABILE STATICA GLOBALE
// ==========================================================
inline Random Random::GLOBAL{Random::STARTUP_SEED};
