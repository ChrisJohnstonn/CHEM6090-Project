function cmass = mass(atomicSymbol)
    switch atomicSymbol
        case 'H'
            cmass = 1.0079;
        case 'He'
            cmass = 4.0026;
        case 'Li'
            cmass = 6.941;
        case 'Be'
            cmass = 9.0122;
        case 'B'
            cmass = 10.811;
        case 'C'
            cmass = 12.0107;
        case 'N'
            cmass = 14.0067;
        case 'O'
            cmass = 15.9994;
        case 'F'
            cmass = 18.9984;
        case 'Ne'
            cmass = 20.1797;
        case 'Na'
            cmass = 22.9897;
        case 'Mg'
            cmass = 24.305;
        case 'Al'
            cmass = 26.9815;
        case 'Si'
            cmass = 28.0855;
        case 'P'
            cmass = 30.9738;
        case 'S'
            cmass = 32.065;
        case 'Cl'
            cmass = 35.453;
        case 'K'
            cmass = 39.0983;
        case 'Ar'
            cmass = 39.948;
        case 'Ca'
            cmass = 40.078;
        case 'Sc'
            cmass = 44.9559;
        case 'Ti'
            cmass = 47.867;
        case 'V'
            cmass = 50.9415;
        case 'Cr'
            cmass = 51.9961;
        case 'Mn'
            cmass = 54.938;
        case 'Fe'
            cmass = 55.845;
        case 'Ni'
            cmass = 58.6934;
        case 'Co'
            cmass = 58.9332;
        case 'Cu'
            cmass = 63.546;
        case 'Zn'
            cmass = 65.39;
        case 'Ga'
            cmass = 69.723;
        case 'Ge'
            cmass = 72.64;
        case 'As'
            cmass = 74.9216;
        case 'Se'
            cmass = 78.96;
        case 'Br'
            cmass = 79.904;
        case 'Kr'
            cmass = 83.8;
        case 'Rb'
            cmass = 85.4678;
        case 'Sr'
            cmass = 87.62;
        case 'Y'
            cmass = 88.9059;
        case 'Zr'
            cmass = 91.224;
        case 'Nb'
            cmass = 92.9064;
        case 'Mo'
            cmass = 95.94;
        case 'Tc'
            cmass = 98;
        case 'Ru'
            cmass = 101.07;
        case 'Rh'
            cmass = 102.9055;
        case 'Pd'
            cmass = 106.42;
        case 'Ag'
            cmass = 107.8682;
        case 'Cd'
            cmass = 112.411;
        case 'In'
            cmass = 114.818;
        case 'Sn'
            cmass = 118.71;
        case 'Sb'
            cmass = 121.76;
        case 'I'
            cmass = 126.9045;
        case 'Te'
            cmass = 127.6;
        case 'Xe'
            cmass = 131.293;
        otherwise
            disp(string('ATOMIC SYMBOL NOT RECOGNISED: ') + atomicSymbol);
    end 