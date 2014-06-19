#!/home/daniel/anaconda/bin/python

def main():
    for i in range(99):
        temp = sing(99-i,'ding dong','alienzzz')
        if(temp == 1):
            break
        if(temp == 2):
            raise NameError('It r tiem 2 dai, todd. Taim 2 dai.\n In other words, the provided input has a value <=0')
    return

def sing(num, item, aliquot):
    if(num > 2):
        print('{} {}s of {} on the wall,\n {} {}s of {}!\n'.format(num,aliquot,item,num,aliquot,item)),
        print('Take one down, pass it around,\n {} {}s of {} on the wall!'.format(num-1,aliquot,item))
        return 0
    elif(num == 2):
        print('2 {}s of {} on the wall,\n 2 {}s of {}!\n'.format(aliquot,item,aliquot, item)),
        print('Take one down, pass it around,\n 1 {} of {} on the wall!'.format(aliquot,item))
    elif(num == 1):
        print('1 {} of {} on the wall,\n just 1 {} of {}!\n'.format(aliquot,item,aliquot,item)),
        print('Take one down, pass it around,\n ... no more {}s of {} on the wall.'.format(aliquot,item))
        return 1;
    else:
        print("What happened? Someone set up us the bomb.")
        return 2;

if __name__ == "__main__":
    main()
