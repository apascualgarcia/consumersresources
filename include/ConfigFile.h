// Classe facilitant la lecture de fichiers de configuration.
// Contributeurs : K. Steiner, J. Dominski, N. Ohana
// Utilisation : Envoyer au constructeur le nom d'un fichier contenant
// les parametres sous la forme [param=valeur] sur chaque ligne, puis
// appeler get<type>("param") pour acceder a un parametre.

#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>


class ConfigFile{

  public:
      ConfigFile(const std::string& filename);
      ~ConfigFile();

      template<typename T> T get(const std::string& key) const;

      void process(const std::string& lineread);

      std::string toString() const;

      void printOut(std::string path) const;

  private:
      std::string trim(std::string str);

      std::map<std::string, std::string> configMap;
};

#endif
