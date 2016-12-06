package utils;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by Liyanzhen on 2016/12/6.
 */
public class RegExpTest {
    public static void main(String [] args){
        String str = "D:\\paperdata\\soybean\\RMA normalization data\\GSE7108_Su_RMA_matrix.txt";

//        Pattern pattern = Pattern.compile("^GSE.*_$");
        MyPrint.print("测试指数求法","" );


        Pattern pattern = Pattern.compile("(GSE){1}[0-9]*");
        Matcher matcher = pattern.matcher(str);
        StringBuffer buffer = new StringBuffer();
        while(matcher.find()){
            buffer.append(matcher.group());
            System.out.println(buffer.toString());
        }
    }
}
