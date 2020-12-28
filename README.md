# TP1 : 

Il y a encore certaine erreur dont le calcul de Rusanov (qui fonctionne pour l'équation du transport mais pas pour burgers).

Cette branche implémente la résolution des équations hyperboliques en 1 dimension 
avec les méthode de Godunov, Rusanov et MUSCL. Elle contient les codes source de deux application :
* solve : résout le problème avec les différentes méthodes
* compute_error : calcul les erreurs entre le résulat des méthodes et la solution exacte en fonction du paramètre N (nombre de points en espace)

## use_me

Un petit fichier bash est présent à la racine qui éxécute les différentes étapes suivantes (compilation puis éxécution de certain problèmes). Aprèes l'éxécution de ce fichier,  les résultats seront présents dans le dossier output.

Attention, son éxécution peut prendre un peu de temps (une dizaine de minute sur mon ordinateur).

## Architecture du projet

Le projet est composé par :

* src : code source
    * solve.c : implemente l'application solve
    * compute_error.c : implemente l'application compute_error
    * parameters.c : contient les classes qui regroupe les différents paramètres
    * parameters_equation.c : contient les différentes fonctions servant aux deux applications
* input : contient des fichiers comprenant les paramètres des problèmes à calculer
* Latex : contient les fichier source du rapport

## Compilation et execution des applications

Nous avons deux applicaiton :
- solve
- compute_error

Pour les compiler vous pouvez utilisez gcc :
- gcc solve.c -o solve -lm
- gcc compute_error.c -o compute_error -lm

Ou cmake en utilisant le fichier CMakeLists.txt :
- créer un fichier build (mkdir build) puis aller dedans
- cmake ..
- make

### Executions

#### solve

Usage de l'application :
    solve {g|r|m|e} path_input path_output

path_input : chemin du fichier d'entré. Ce fichier contient les paramètres du problème :
    - option_godunov : equation à résoudre
    - xmin
    - xmax
    - cfl
    - N : nombre de point en espace
    - tmax : borne temporelle
    Vous pouvez voir des examples dans input

path_output : chemin du dossier ou sera créer le dossier de sortie (du même nom que le fichier d'entrée) qui contiendra les résultats. Il est structuré comme ceci :
    - parameters : contient les différents paramètres du problème
    - plot.dat : résultat du problème en espace
    - plotcom.gnu : sert à générer le graphique
    - plot.png : graphique des résultats en espace

option :
    - g résout le problème avec la méthode de Godunov
    - r résout le problème avec la méthode de Rusanov
    - m résout le problème avec la méthode MUSCL
    - e calcul l'erreur exacte

#### compute_error

Usage de l'application :
    compute_error {g|r|m} path_input path_output

path_input : chemin du fichier d'entré. Ce fichier contient les paramètres du problème :
    - option_error : norme de l'erreur
    - option_godunov : equation à résoudre
    - xmin
    - xmax
    - cfl
    - len_liste_N : longueur de liste_N
    - liste_N : liste des dimension spatiale N servant à calculer l'erreur
    - tmax : borne temporelle
    Vous pouvez voir des examples dans input

path_output : chemin du dossier ou sera créer le dossier de sortie qui contiendra les résultats. Il est structuré comme ceci :
    - parameters : contient les différents paramètres du problème
    - plot.dat : résultat du calcul (temps et erreur)
    - plotcom.gnu : sert à générer le graphique
    - error.png : graphique des erreur en fonction de N
    - time.png : graphique des durée en fonction de N

option :
    - g calcul l'erreur avec la méthode de Godunov
    - r calcul l'erreur avec la méthode de Rusanov
    - m rcalcul l'erreur avec la méthode MUSCL

