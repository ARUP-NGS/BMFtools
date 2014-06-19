#/home/daniel/anaconda/bin/python

def main():
    for i in range(99):
        temp = sing(99-i)
        if(temp == 1):
            break
        if(temp == 2):
            raise NameError('It r tiem 2 dai, todd. Taim 2 dai.\n In other words, the provided input has a value <=0')
    return

def sing(num):
    if(num > 1):
        print('{} bottles of beer on the wall,\n {} bottles of beer!\n'.format(num,num)),
        print('Take one down, pass it around,\n {} bottles of beer on the wall!'.format(num-1))
        return 0
    elif(num == 1):
        print('1 bottle of beer on the wall,\n just 1 bottle of beer!\n'),
        print('Take one down, pass it around,\n ... no more bottles of beer on the wall.')
        return 1
    else:
        print("What happened? Someone set up us the bomb.")
        return 2

if __name__ == "__main__":
    main()
