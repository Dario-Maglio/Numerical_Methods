#ifndef PHONEBOOK_H
#define PHONEBOOK_H
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <utility>
#include <fstream>
#include <sstream>

class PhoneBookIN {
public:
  PhoneBookIN() {}
  std::string name, telnum, email;
  std::map<std::string, std::string> addinfo;
  PhoneBookIN(const std::string &nam, const std::string &tel,
                                  const std::string &mail=""):
    name(nam), telnum(tel), email(mail) {}

  void print() const {
    std::cout<<name<<" "<<telnum<<" "<<email<<std::endl;
    for (std::map<std::string,std::string>::const_iterator it = addinfo.begin();
                                                    it != addinfo.end(); it++) {
      std::cout<< " -"<< it->first << " is "<< it->second<<std::endl;
    }
  }

  bool operator < (const PhoneBookIN &other) const {
    return (name < other.name);
  }
};

//------------------------------------------------------------------------------

class PhoneBook {
private:
  bool sw_;

public:
  PhoneBook():
    sw_(false) {}
  std::vector<PhoneBookIN> contacts;

  void addCont(const std::string &nam, const std::string &tel,
               const std::string &mail ="") {
    PhoneBookIN contact(nam, tel, mail);
    contacts.push_back(contact);
    sw_ = false;
  }

  void addPhoneCont(const PhoneBookIN &phoneIn) {
    contacts.push_back(phoneIn);
    sw_ = false;
  }

  void printContacts(){
      for (size_t i=0; i<contacts.size(); i++) {
        contacts[i].print();
      }
  }

  /***SORTING & FINDING METHODS***/

  void sorter(){
    std::sort(contacts.begin(), contacts.end());
    sw_ = true;
  }

  int recurfinder(const std::string &name, const int &a, const int &b) {
    int c = (b + a)/2;
    if (name == contacts[c].name) {
      return c;
    } else if (name < contacts[c].name) {
      return recurfinder(name, a, c);
    } else if (c != a) {
      return recurfinder(name, c, b);
    }
    //Se arriva qui, allora è in loop tra due elementi
    if (name == contacts[b].name) {
      return b;
    } else {
      return -1;
    }
  }

  std::pair<bool, int> finder(const std::string &name) {
    if (sw_ == false){
      sorter();
    }
    int c = recurfinder(name, 0, contacts.size() - 1);
    if (c != -1) {
      std::cout<<"Element "<<name<<" found in position "<< c <<std::endl;
      return std::make_pair(true, c);
    } else {
      std::cout<<"Element "<<name<<" not in contacts!"<<std::endl;
      return std::make_pair(false, c);
    }
  }

  /***TRANSFORMING METHODS***/

  static PhoneBookIN prefAdd(const PhoneBookIN &contact) {
    PhoneBookIN newcontact = contact;
    if (contact.telnum == "" or contact.telnum[0] != '+') {
      newcontact.telnum = "+39" + contact.telnum;
    }
    return newcontact;
  }

  void addPrefix () {
    std::transform(contacts.begin(), contacts.end(), contacts.begin(), prefAdd);
  }

  /***FILES METHODS***/

  void loadPhoneBook(const std::string &path) {
    //Load a file with format "name, telephon, email"
    std::ifstream file(path);
    if (file.is_open()) {
      std::string line;
      while (getline(file,line)){
        std::string field;
        std::vector<std::string> fields;
        std::stringstream strmstr(line);
        while (getline(strmstr, field, ',')) fields.emplace_back(field);
          addCont(fields[0], fields[1], fields[2]);
      }
      file.close();
    } else std::cout << "Unable to open file."<<std::endl;
  }

  void savePhoneBook(){
    //Save a file with format "name, telephon, email"
    int len = contacts.size();
    if (len != 0) {
      std::ofstream file;
      file.open ("Phone book.txt");
      for (size_t i=0; i<len; i++) {
        file <<contacts[i].name<<", "<<contacts[i].telnum<<", "<<
                                    contacts[i].email<<std::endl;}
      file.close();
    }
  }
};

#endif
