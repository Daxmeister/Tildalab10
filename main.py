# Lab 10 Vira Oetterli och Davide Attebrant Sbrzesny

from linkedQFile import LinkedQ
from molgrafik import *
from hashtable_atom import *
#######################################################################################################################
# Stödfunktioner
#######################################################################################################################


def create_list_of_atoms():
    '''Returnerar en lista med alla kända atomer.'''
    string = 'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr ' \
             'Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In ' \
             'Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re ' \
             'Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md ' \
             'No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Fl Lv'
    list_of_atoms = string.split(' ')
    return list_of_atoms


def enqueue_formel(formel_string):
    '''Konverterar en sträng till en linked queue.'''
    queue = LinkedQ()
    for character in formel_string:
        if character != '\n' and character != ' ':
            queue.enqueue(character)
    return queue

class Ruta:
    '''Skapar en ruta som är till för en grupp eller parentes'''
    def __init__(self, atom="( )", num=1):
        self.atom = atom
        self.num = num
        self.next = None
        self.down = None

    def return_weight(self):
        atom_dict = create_atom_dict()
        if self.atom == "( )":
            print(self.num)
            return
        else:
            return atom_dict[self.atom] * int(self.num)



#######################################################################################################################
# BNF-Syntax functioner
#######################################################################################################################


class Syntaxfel(Exception):
    pass


def read_formel(queue):
    '''<formel>::= <mol> \n'''
    return read_molekyl(queue)


def read_molekyl(queue):
    '''<mol>   ::= <group> | <group><mol>'''

    mol = read_group(queue) # Skapar ett "ruta" klass objekt
    if not queue.isEmpty():     # Den fortsätter läsa in grupper tills den stöter på en slutparentes.
        if queue.peek() == ")":
            global number_open_paranthesis_global   # Denna variabel håller koll på om vi har öppnat några parenteser.
            if number_open_paranthesis_global == 0: # Har vi inga öppnade parenteser ska vi inte få slutparenteser.
                raise Syntaxfel("Felaktig gruppstart vid radslutet ")
            else:
                return mol # Detta tar oss tillbaka till där vi kallades ifrån. Antingen kolla_molekyl eller
                        # till read_group efter att en parentes har öppnats. Där kommer nästa steg att kolla siffror.
        mol.next = read_molekyl(queue)
    return mol

# Grupp. Får börja med atom eller start-parentes.
# Är gruppstart inte en bokstav eller "(" ska vi få meddelande fel gruppstart.
# Annars ska vi antingen läsa in atom + nummer
# Eller ropa på molekyl igen och sen läsa nummer.
def read_group(queue):
    '''<group> ::= <atom> |<atom><num> | (<mol>) <num>'''
    # Skapa ruta
    rutan = Ruta()


    # Kollar först att vi startar gruppen med godkänd karaktär.
    first_character = queue.peek()
    acceptable_groupstart_characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz("
    if first_character not in acceptable_groupstart_characters:
        raise Syntaxfel("Felaktig gruppstart vid radslutet ")

    # Om det inte är en (<mol>) start läser den in den som en atom.
    if first_character != "(":
        '''<atom> | <atom><num>'''
        rutan.atom = read_atom(queue)
        numbers = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
        if queue.peek() in numbers:
            rutan.num = read_number(queue)
        return rutan

    # Om vi börjar med "(" kommer den att hantera det som en molekyl.
    else:
        '''(< mol >) < num >'''
        rutan.atom = "( )"
        queue.dequeue() # Ta bort startparentes
        global number_open_paranthesis_global
        number_open_paranthesis_global += 1 # Lägger till en öppnad parentes till vår global variabel.
        rutan.down = read_molekyl(queue)
        if queue.peek() != ")":
            raise Syntaxfel("Saknad högerparentes vid radslutet ")
        queue.dequeue() # Tar bort slutparentes
        number_open_paranthesis_global -= 1 # Tar bort en öppnad parentes till vår global variabel.

        # Kollar att det kommer en siffra efter avslutad parentes
        numbers = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
        if queue.peek() not in numbers:
            raise Syntaxfel ("Saknad siffra vid radslutet ")
        rutan.num = read_number(queue)
        return rutan



# Kollar om stor bokstav, ger fel annars. Kollar om nästa är en liten bokstav och kommer isåfall tugga upp den också.
# Jämför dessa två med bekanta atomer. OBS att den inte kollar vidare på queue. Läser max in två tecken.
# Klarar inte att få None skickat till sig.
def read_atom(queue):
    '''<atom>  ::= <LETTER> | <LETTER><letter> OCH returnerar atomen som string'''
    lowercase_letters = "abcdefghijklmnopqrstuvwxyz"
    capital_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    # Kollar att första bokstaven är stor, annars fel
    next_character = queue.peek()
    if next_character not in capital_letters:
        raise Syntaxfel("Saknad stor bokstav vid radslutet ")

    #Läser första bokstaven
    read_capital_letter(queue)

    # Kollar på nästa tecken, är det en liten bokstav blir den inläst.
    next_next_character = queue.peek()
    atom = next_character
    if next_next_character != None and next_next_character in lowercase_letters:
        atom = next_character + next_next_character
        read_lowercase_letter(queue)

    # Jämför vår atom av en eller två bokstäver med listan.
    accepted_atoms_list = create_list_of_atoms()
    if atom in accepted_atoms_list:
        return atom
    else:
        raise Syntaxfel("Okänd atom vid radslutet ")



# Kan få vad som helst, ser om stor bokstav
def read_capital_letter(queue):
    '''<LETTER>::= A | B | C | ... | Z'''
    capital_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    character = queue.dequeue()
    if character in capital_letters:
        return
    else:
        raise Syntaxfel("Saknad stor bokstav vid radslutet " + character)


# Kan få vad som helst, ser om det är liten bokstav.
def read_lowercase_letter(queue):
    '''<letter>::= a | b | c | ... | z'''
    lowercase_letters = "abcdefghijklmnopqrstuvwxyz"
    character = queue.dequeue()
    if character in lowercase_letters:
        return
    else:
        raise Syntaxfel("Saknad stor bokstav vid radslutet ")


def read_number(queue):
    '''<num>   ::= 2 | 3 | 4 | ... OCH returnerar siffran som int'''
    number_string = queue.dequeue()
    if number_string == "0":
        raise Syntaxfel("För litet tal vid radslutet ")
    numbers = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
    while queue.peek() in numbers:  # Alla siffror på rad dequeueas.
        number_string = number_string + queue.dequeue()
    number_int = int(number_string)
    if number_int < 2:
        raise Syntaxfel("För litet tal vid radslutet ")
    return number_int

'''test_string = "1"
test_queue = enqueue_formel(test_string)
print(read_number(test_queue))
exit()'''


#######################################################################################################################
# Molekylvikt
#######################################################################################################################
def weight(first_ruta): # Horizontal approach
    '''Använder sig av trödet med rutor för att räkna ut vikt.'''
    return horizontal_weight(first_ruta)


def horizontal_weight(ruta):
    # Kollar om det finns ne till åt "höger" och går isåfall dit. Annars kollar den om den kan ta vikten rakt av och
    # addera eller om den ör vid en parentes. Är det en parentes måste den jobba sig ned vertikalt.
    if ruta.next == None:   # Vi är vid slutstationen.
        if ruta.atom == "( )":
            print("A")
            print(ruta.num)
            return horizontal_weight(ruta.down) * int(ruta.num)

        else:
            print("B")
            return ruta.return_weight()

    else:   # Vi är vid en mellanstation.
        if ruta.down == None:
            print("C")
            return horizontal_weight(ruta.next) + ruta.return_weight()

        # Detta betyder att vi är en parentes
        else:
            print("D")
            return horizontal_weight(ruta.next) + horizontal_weight(ruta.down) * int(ruta.num)



def weight_ver(ruta):
    return vertical_weight(ruta)

def vertical_weight(ruta):
    if ruta.down == None: # Om det inte finns något nedanför. Dvs. ingen parentes.
        return weight(ruta)
    else:  # Om det finns något nedanför. Dvs. vi har en parentes.
        return vertical_weight(ruta.down) * int(ruta.num)



#######################################################################################################################
# Run funktioner
#######################################################################################################################


# Input är en string, output är ett meddelande om ifall den följer syntax eller inte.
def kolla_molekyl(molekylstring):
    '''Tar emot en string och kör programmet.'''
    queue = enqueue_formel(molekylstring)
    global number_open_paranthesis_global
    number_open_paranthesis_global = 0
    try:
        mg = Molgrafik()
        first_ruta = read_formel(queue)
        mg.show(first_ruta)  # Här ritar vi ut den.
        print(weight(first_ruta))
        input("redo för nästa")
        return "Formeln är syntaktiskt korrekt"
    except Syntaxfel as fel:
        return str(fel) + str(queue)


def main():
    '''Kör själva programmet. Ger möjlighet för input och skickar till kolla_molekyl.'''
    string_molekyl = input()
    while string_molekyl != '#':
        print(kolla_molekyl(string_molekyl))
        string_molekyl = input()

main()