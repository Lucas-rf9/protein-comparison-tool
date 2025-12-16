import argparse
#lib qui permet d'utiliser le code en ligne de commandes


#Partie 1: lecture du fichier fasta
def extraire_prot_fasta(fichier_fasta):
    #on lit un fichier fasta ligne par ligne et on récupère les séquences de protéines
    proteines = {}
    a = None  #variable qui sert à stocker l'identifiant actuel

    for ligne in fichier_fasta:
        ligne = ligne.strip('\n')  #on enlève le retour
        if not ligne:
            continue  #on saute les lignes vides dans le cas où il y en aurait

        if ligne[0] == '>':
            #on s'intéresse aux lignes qui commencent par '>' car elles contiennent l'identifiant de la protéine
            a = ligne[4:10]  #on récupère l'ID de la prot
            proteines[a] = ""  #on initialise une clé dans le dict pour cette prot
        else:
            #sinon c'est une ligne de séquence, on rajoute les aa à la suite
            proteines[a] += ligne

    return proteines


#Partie 2: création du dictionnaire de k-mers
def seq_to_kmers(sequence_prot, k):
    #la fonction prend une séquence protéique et un entier k comme arguments
    #elle retourne un dict contenant tous les k-mers trouvés + leur nb d'occurrences
    dict_kmers = {}

    #on récupère dans une var une seq de taille k
    for i in range(0, len(sequence_prot) - k + 1):
        kmers = sequence_prot[i : i + k] 
        #si la seq existe déjà on rajoute +1 à sa valeur dans le dict
        if kmers in dict_kmers:
            dict_kmers[kmers] += 1
        #sinon on initialise la valeur de sa clé à 1
        else:
            dict_kmers[kmers] = 1

    return dict_kmers


#Partie 3: calcul de la proportion de k-mers uniques
def proportion_kmers_uniques(dict_seq_A, dict_seq_B):
    #on convertit les clés, kmers, en ensembles cad on retire les doublons de clés et on regarde si ils ont des clés en commun
    ensemble_A = set(dict_seq_A.keys())
    ensemble_B = set(dict_seq_B.keys())

    #nombre de kmers uniques en commun
    kmers_commun = len(ensemble_A & ensemble_B)

    #calcul du dénominateur = le plus petit nb de kmers uniques entre les deux séquences
    denominateur = min(len(ensemble_A), len(ensemble_B))

    #calcul de la proportion
    return kmers_commun / denominateur


#Partie 4: calcul de la proportion de kmers totaux 
def proportion_kmers_commun(dict_A, dict_B):
    #on veut savoir quelle est la proportion de kmers en commun de deux seq
    #en prenant en compte le nb total, avec les doublons
    
    communs = 0  #la var stocke combien de kmers sont communs au total
    
    #on parcourt tous les kmers de la première seq
    for kmer in dict_A:
        #si le kmer existe aussi dans la 2e seq
        if kmer in dict_B:
            #on ajoute le plus petit nb d’occurrences entre les deux
            communs += min(dict_A[kmer], dict_B[kmer])

    #on récupère le total de kmers pour chaque séquence en faisant la somme des valeurs de chaque dict
    totalA = sum(dict_A.values())
    totalB = sum(dict_B.values())

    #on calcule la proportion par rapport au plus petit total
    proportion = communs / min(totalA, totalB)

    return proportion


#Partie 5: comparaison de la seq requête avec toutes les seq de la banque
def comparaison_sequences(proteines, id_requete, k):
    #on transforme la seq requête en dict de kmers
    kmers_requete = seq_to_kmers(proteines[id_requete], k)
    resultats = []  #on stocke ici les infos de comparaison pour chaque cible

    #on parcourt toutes les séquences dans la banque
    for id_cible in proteines:
        #dans le cas où la requete est aussi dans la bq de seq on la passe
        if id_cible == id_requete:
            continue

        #on génère les kmers de la séquence cible
        kmers_cible = seq_to_kmers(proteines[id_cible], k)

        #comparaison sur les kmers totaux
        communs = 0  #compteur de k-mers en commun 
        for kmer in kmers_requete:
            if kmer in kmers_cible:
                communs += min(kmers_requete[kmer], kmers_cible[kmer])

        #total de kmers dans chaque séquence
        total_requete = sum(kmers_requete.values())
        total_cible = sum(kmers_cible.values())

        #on divise le nb de k-mers uniques en commun par le plus petit total de kmers uniques des deux
        prop_totaux = communs / min(total_requete, total_cible)

        #comparaison sur les kmers uniques
        set_requete = set(kmers_requete.keys())
        set_cible = set(kmers_cible.keys())

        #on compte combien de kmers uniques sont partagés
        uniques = len(set_requete & set_cible)
        
        prop_uniques = uniques / min(len(set_requete), len(set_cible))

        #on stocke les infos dans un dict pour cette séquence cible
        resultats.append({
            "id": id_cible,
            "nb_kmers_totaux": communs,
            "prop_kmers_totaux": prop_totaux,
            "nb_kmers_uniques": uniques,
            "prop_kmers_uniques": prop_uniques
        })

    #tri des résultats
    #on trie pour avoir les séquences les plus similaires en premier
    def tri_par_proportion(dict):
        return dict["prop_kmers_totaux"]
    
    #valeur de tri utilisée : proportion de kmers totaux
    resultats.sort(key=tri_par_proportion, reverse=True)

    #on renvoie la liste finale triée
    return resultats


#Partie 6: écriture du fichier de résultats
def ecriture_fichier_compte(liste_resultats, nom_fichier):
    #on crée un fichier texte tabulé avec les résultats de la comparaison
    with open(nom_fichier, "w") as f:
        #header du tableau
        f.write("id\tnb_kmers_totaux\tprop_kmers_totaux\tnb_kmers_uniques\tprop_kmers_uniques\n")

        #on écrit les infos pour chaque séquence cible
        for dico in liste_resultats:
            f.write(dico["id"] + "\t" +
                    str(dico["nb_kmers_totaux"]) + "\t" +
                    str(dico["prop_kmers_totaux"]) + "\t" +
                    str(dico["nb_kmers_uniques"]) + "\t" +
                    str(dico["prop_kmers_uniques"]) + "\n")


#Partie 7: écriture des kmers de la requête
def ecriture_kmers_requete(dict_kmers_requete, nom_fichier):
    #on veut enregistrer tous les kmers trouvés dans la requête
    #en utilisant leur nombre d’occurrences, triés du plus fréquent au moins fréquent

    with open(nom_fichier, "w") as f:
        f.write("kmer\tcompte\n")

        #on doit convertir le dict en liste de tuples pour pouvoir le trier
        liste_kmers = list(dict_kmers_requete.items())

        #on tri par fréquence décroissante
        def trier_par_freq(element):
            return element[1]  #on renvoie la fréquence 
        
        liste_kmers.sort(key=trier_par_freq, reverse=True)

        #ecriture de chaque kmer dans le fichier
        for kmer, nb in liste_kmers:
            f.write(kmer + "\t" + str(nb) + "\n")


#Partie 8: affichage d’un petit résumé en console
def afficher_resume(id_requete, dict_kmers_requete, liste_resultats):
    #on affiche combien de k-mers uniques contient la requête
    nb_uniques = len(set(dict_kmers_requete.keys()))
    print(f"{id_requete} contient {nb_uniques} k-mers uniques.")

    #on affiche la séquence la plus proche
    if len(liste_resultats) > 0:
        meilleure = liste_resultats[0]
        print("La séquence la plus proche possède " + str(meilleure["nb_kmers_totaux"]) +
              " k-mers en commun, avec l’identifiant " + meilleure["id"] + ".")
    else:
        print("Aucune séquence cible trouvée.")


#Fonction principale : Main fonction
if __name__ == "__main__":
    #argparse permet de lancer le script dans le terminal
    parser = argparse.ArgumentParser()
    #ajout des 3 arguments de la fonction: la séquence requête, la banque de séquences et la taille des kmers
    parser.add_argument("fichier_requete")  
    parser.add_argument("fichier_banque")   
    parser.add_argument("k", type=int)   
    args = parser.parse_args()

    #lecture du fichier contenat la banque de seq
    with open(args.fichier_banque, "r") as f_banque:
        proteins = extraire_prot_fasta(f_banque.readlines())

    #lecture du fichier de la seq requête
    with open(args.fichier_requete, "r") as f_req:
        dict_requete = extraire_prot_fasta(f_req.readlines())

    #on récupère l'identifiant de la seq requete
    id_requete = list(dict_requete.keys())[0]

    #on ajoute la seq requete dans la banque pour pouvoir la comparer
    proteins[id_requete] = dict_requete[id_requete]

    #création du dictionnaire de kmers pour la requête
    dict_kmers_requete = seq_to_kmers(proteins[id_requete], args.k)

    #on compare la séquence requête avec toutes les séquences cibles
    comp_seq = comparaison_sequences(proteins, id_requete, args.k)

    #écriture des résultats dans les fichiers de sortie
    ecriture_kmers_requete(dict_kmers_requete, "kmers_requete.tsv")
    ecriture_fichier_compte(comp_seq, "comparaison_resultats.tsv")

    #affichage de l'output dans le terminal
    afficher_resume(id_requete, dict_kmers_requete, comp_seq)
