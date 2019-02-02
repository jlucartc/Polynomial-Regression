from sympy import * ;

import time;

import numpy as np;

from math import * ;

# PT-BR : monta a equação de ordem fornecida com os símbolos para posterior substituição
# EN : builds up the equation of n-th order with symbols for further substitution

def montarEquacao(incognitas,x): 
    
    y = 0;

    for i in range(0,len(incognitas)):

        y += incognitas[i]*(x**(len(incognitas)-(i+1)));
        
    return y;


# PT-BR: checa se o valor fornecido pode ser interpretado como float
# EN : checks if the given value can be seen as float

def isFloat(s) : 

    try :

        float(s);
        return True;

    except:

        return False;

'''
PT-BR : A ordem máxima vai até o 24º grau pois para 24 graus são necessárias 25
        incognitas, sendo que o alfabeto possui 26 letras e 1 está reservada p/ substituição (a letra 'x')
'''

'''
EN : The limit goes until the 24-th order, because the alphabet only has 26
     letters, from which one is reserved for value substitution (the letter 'x'),
     and the remaining 25 can represent at most a 24-th order polynomial.
'''

n = int(input("Digite a ordem da regressão (ordem máxima = 24): ")); 

if(n > 24):

    print("\n\n------------------\n  Ordem inválida  \n------------------\n\n");
    
else:

    #PT-BR: variável reservada p/ substituição
    #EN: reserved variable for substitution

    x = symbols('x');

    # PT-BR : abre o arquivo csv onde os dados se encontram
    # EN: opens the csv file where the data is saved
    
    f = open('data.csv','r');

    #PT-BR: lista de simbolos disponíveis para representar os coeficientes da curva
    #EN : list of symbols available for the curve's coefficients representation

    lista_simbolos = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','y','z'];

    # PT-BR: lista das abcissas
    # EN: abcissas' list

    xLista = [];

    # PT-BR: lista das ordenadas
    # EN: ordenates' list
    
    yLista = []; 


    # PT-BR: insere os valores nos vetores x e y
    # EN: inserts the values in x and y arrays

    for i in f: 

        s = [];
        s = i.split(',',1);
        s[1] = s[1].split('\n',1)[0];
        
        if(isFloat(s[0]) and isFloat(s[1])):

            xLista.append(float(s[0]));
            yLista.append(float(s[1]));

        else:

            continue;

    incognitas = [];

    # PT-BR: cria o vetor com as incognitas que serão utilizadas
    # EN: creates the vector with all the incognites that will be used

    for i in range(0,n+1): 

        incognitas.append(symbols(lista_simbolos[i]));

    
    h = montarEquacao(incognitas,x);

    somaQuadradaDiferencas = 0;

    # PT-BR: realiza a soma dos quadrados dos erros em relação a cada ponto
    # EN: makes the summation of the squared errors with respect to each given point
    
    for i in range(0,len(xLista)): 

        somaQuadradaDiferencas += ( yLista[i] - h.subs(x,xLista[i]))**2


    diffs = [];

    # PR-BR: calcula a derivada da soma dos quadrados com relação a cada incognita
    # EN : takes the derivative of the summation of squared errors relative to each incognite

    for i in range(0,len(incognitas)): 

        diffs.append(diff(somaQuadradaDiferencas,incognitas[i]));


    rep = [];

    # PT-BR: resolve cada derivada em relação a cada incognita
    # EN: solves the derivatives by isolating each incognite
    
    for i in range(0,len(incognitas)): 
    
        if i == 0 :

            diffs[len(incognitas)-1] = solve(diffs[len(incognitas)-1],incognitas[len(incognitas)-1])[0];

            rep.append((incognitas[len(incognitas)-1],diffs[len(incognitas)-1]));
            
        else :

            diffs[len(incognitas)-(i+1)] = solve(diffs[len(incognitas)-(i+1)].subs(rep),incognitas[len(incognitas)-(i+1)])[0];

            rep.append((incognitas[len(incognitas)-(i+1)],diffs[len(incognitas)-(i+1)]));


    rep = [(incognitas[0],diffs[0])];


    # PT-BR: substitui os resultados nas equações subsequentes p/ obter o valor de cada incognita
    # EN: makes the substitution of the results on the subsequent equations so the icognites' values can be found
    
    for i in range(0,n): 

        diffs[i+1] = diffs[i+1].subs(rep);
        rep.append((incognitas[i+1],diffs[i+1]));

    

    k = 0;

    # PT-BR: monta a expressão final
    # EN: builds the final expression
    
    for i in range(0,len(incognitas)):
        
        k += diffs[i]*x**(len(incognitas)-(i+1))

    print("Modelo: ",k);
