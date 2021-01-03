# TP2 : 

Il y de gros problème pour la résolution de Saint Venant (pour la méthode de Godunov et de Rusanov).

Cette branche implémente la résolution de l'équation de Saint Venant en dimension 1 avec les méthode de Godunov, Rusanov et MUSCL. Elle contient les codes source de deux applications :
* solve : résout le problème avec les différentes méthodes
* riemann : calcul le solveur de Riemann et le plot (elle n'est pas compilé par me cmake).

## use_me

Un petit fichier bash est présent à la racine qui éxécute les différentes étapes suivantes (compilation puis éxécution de certain problèmes). Aprèes l'éxécution de ce fichier,  les résultats seront présents dans le dossier output.

## Architecture du projet

Le projet est composé par :

* src : code source
    * riemann : contient les code de l'application riemann
    * solve : contient les codes de l'application solve
* input : contient des fichiers comprenant les paramètres des problèmes à calculer
* Latex : contient les fichier source du rapport

## Compilation et execution des applications

Nous avons deux applicaiton :
- solve
- riemann

Pour les compiler vous pouvez utilisez gcc :
- gcc solve.c -o solve -lm
- gcc riemann.c -o riemann -lm

Ou cmake en utilisant le fichier CMakeLists.txt :
- créer un fichier build (mkdir build) puis aller dedans
- cmake ..
- make

Remarque : riemann n'est pas compilé par cmake mais solve oui.

### Executions

#### solve

Usage de l'application :
    solve -gr path_input path_output

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
