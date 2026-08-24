#pragma once

#include <Rinternals.h>

// Includi il file core se necessario in questo header,
// altrimenti puoi includerlo direttamente nel file .cpp
#include "separation.h"

/**
 * @brief Functor per eseguire funzioni R personalizzate (t-norm e t-conorm).
 * * Questo functor fa da ponte tra il loop HPC in puro C++ e l'interprete R.
 * Utilizza la pre-allocazione dei SEXP (arg1, arg2) per evitare le allocazioni
 * di memoria durante l'esecuzione del ciclo intensivo, abbattendo l'overhead del
 * Garbage Collector di R.
 * * Soddisfa il concept `NormConormFunc` definito in separation.h.
 */
struct RCustomNormConorm {
    SEXP arg1;
    SEXP arg2;
    SEXP call_expr;
    SEXP env;
    
    // Costruttore
    RCustomNormConorm(SEXP a1, SEXP a2, SEXP call, SEXP e)
    : arg1(a1), arg2(a2), call_expr(call), env(e) {}
    
    // Overload dell'operatore () per simulare una funzione C++
    inline double operator()(double a, double b) const {
        // 1. Modifica i valori in-place nei vettori R pre-allocati
        REAL(arg1)[0] = a;
        REAL(arg2)[0] = b;
        
        // 2. Valuta l'espressione (la funzione R) nell'ambiente fornito
        SEXP res = Rf_eval(call_expr, env);
        
        // 3. Estrae e restituisce il double
        return Rf_asReal(res);
    }
};
