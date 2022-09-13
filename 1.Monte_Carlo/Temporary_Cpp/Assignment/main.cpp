#include <iostream>
#include "phonebook.h"

int main() {
  PhoneBook rubrica;
  rubrica.addCont("Giovanni", "3286243153", "abc@unipi.it");
  rubrica.addCont("Miovi", "");
  PhoneBookIN contact("Aiovi", "+39329", "hafg@pg");
  contact.addinfo["Class"] = "2L";
  contact.addinfo["Address"] = "AV";
  PhoneBookIN contact2("Bruno", "3286243153", "hag@gma");
  contact2.addinfo["Class"] = "3B";
  rubrica.addPhoneCont(contact);
  rubrica.addPhoneCont(contact2);
  rubrica.printContacts();

  rubrica.savePhoneBook();

  rubrica.finder("Ciccio");
  rubrica.finder("Miovi");
  rubrica.finder("Aiovi");
  rubrica.addPrefix();
  rubrica.printContacts();

  PhoneBook pagine;
  pagine.loadPhoneBook("Phone book.txt");
  pagine.printContacts();
}
