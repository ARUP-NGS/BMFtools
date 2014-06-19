#!/home/daniel/anaconda/bin/python

def main():
    for i in range(99):
        count = 99-i
        print('{} bottles of beer on the wall,\n {} bottles of beer!\n'.format(count,count)),
        print('Take one down, pass it around,\n {} bottles of beer on the wall!'.format(count-1))

if __name__ == "__main__":
    main()
