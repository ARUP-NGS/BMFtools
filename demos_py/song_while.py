#!/home/daniel/anaconda/bin/python

def main():
        count=99 #This is a comment

        while(count > 0):
            print('{} bottles of beer on the wall,\n {} bottles of beer!\n'.format(count,count)),
            count-=1
            print('Take one down, pass it around,\n {} bottles of beer on the wall!'.format(count))

if __name__ == "__main__":
    main()
