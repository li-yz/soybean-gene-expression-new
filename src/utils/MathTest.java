package utils;

import java.util.*;

/**
 * Created by Liyanzhen on 2016/12/6.
 */
public class MathTest {
    public static void main(String [] args){
        double r = Math.pow(2,3);
        MyPrint.print("测试指数求法","" + r);

        double [] d1 = {2.5,2.55,5,5,4,4};
        double[] d2 = {1.2,1.3,1.3,1.32,1.4,1.42};
        double sum1=0.0d;
        double sum2=0.0d;
        double mm=0.0d;
        for(int i=0;i < 6;i=i+2){
            sum1 +=d1[i] * d1[i];
            sum2 +=d2[i] * d2[i];
            mm +=d1[i] *d2[i];
        }

        double cos1 = mm/(Math.sqrt(sum1) * Math.sqrt(sum2));

        sum1=0.0d;
        sum2=0.0d;
        mm=0.0d;
        for(int i=0;i < 6;i=i+1){
            sum1 +=d1[i] * d1[i];
            sum2 +=d2[i] * d2[i];
            mm +=d1[i] *d2[i];
        }

        double cos2 = mm/(Math.sqrt(sum1) * Math.sqrt(sum2));

        MyPrint.print("1:"+cos1,"\t 2 :"+cos2);
    }
}
