# Partial derivative partial
TP of classes EDP2 of Master CSMI


## Introduction

Le projet est l'implementation du schéma de godunov appliqué aux l'équations de transport et de burgers. Il suit la trame donné par le sujet ph-tp1.pdf.

## Architecture du projet

Le projet est composé par :

- src : code source
-- raler.c : implemente une fonction de message d'erreur
-- godunov.h : header des code source
-- parameters.c : implemente les différentes fonctions qui paramètre le problème (solution exacte, norme L1 ...)
-- function.c : implemente les différentes fonctions d'initialisation, de plot,de résolution du schéma de godunov et du calcul de l'erreur en fonction du nombre N
-- godunov.c : implemente le calcul du schéma de godunov
-- godunov_error.c : implemente le calcul de l'erruer entre la résultat du schéma de godunov et la solution exacte
- input : contient des fichiers comprenant les paramètres des problèmes calculés
- output : contient les résultats des résolutions
- Latex : contient les fichier source du rapport

## Compilation et execution des applications

Nous avons deux applicaiton :
- godunov
- godunov_error

Pour les compiler vous pouvez utilisez gcc :
- gcc godunov.c -o godunov -lm
- gcc godunov_error.c -o godunov_error -lm

Ou cmake en utilisant le fichier CMakeLists.txt.

### Execution godunov

Pour executer l'application godunov, nous avons la synthaxe :
- ./godunov init_file output_path

Avec init_file un fichier composé des paramètres du problème. Vous pouvez voire les exmaples dans input.

Avec output_path le chemin du repertoire du dossier contenant les résultats de la résolution. Le nom de ce dossier est le nom du fichier init_file. Ce dossier est composé de la manière suivante :
- parameters : récapitulant touts les paramètres du problèmes
- plot.dat : composés de trois colonnes -> les abscisses, la solution numérique et la solution exacte
- plotcom.gnu : permettant de tracer le graphique des données de plot.dat
- graphe.png : le graphe des données de plot.dat

Ce dossier, ainsi que ses fichiers, sont créés directement par l'application. Le plotcom est vide, il faut le remplir manuellement (puis soit vous réexécuter l'application soit vous tracer directement le graphe en ligne de commande).

### Execution godunov_error

Pour executer l'application godunov, nous avons la synthaxe :
- ./godunov_error init_file output_path

Avec init_file un fichier composé des paramètres du problème. Vous pouvez voire les examples dans input.

Avec output_path le chemin du repertoire du dossier contenant les résultats de la résolution. Le nom de ce dossier est le nom du fichier init_file. Ce dossier est composé de la manière suivante :
- parameters : récapitulant touts les paramètres du problèmes
- error.dat : composés de trois colonnes -> les abscisses, l'erreur
- error.png : le graphe des données de error.dat
- time.dat : composés de trois colonnes -> les abscisses, la durée des résolutions
- time.png : le graphe des données de time.dat
- plotcom.gnu : permettant de tracer le graphique des données de plot.dat

Ce dossier, ainsi que ses fichiers, sont créés directement par l'application. Le plotcom est vide, il faut le remplir manuellement (puis soit vous réexécuter l'application soit vous tracer directement le graphe en ligne de commande).

L'application initialise et résout len_N fois le même problème mais en faisant varier le nombre N (nombre d'abscisses).
